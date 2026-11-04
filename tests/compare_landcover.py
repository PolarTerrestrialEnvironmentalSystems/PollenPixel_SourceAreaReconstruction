"""Compare final land-cover stacks against another reconstruction.

Requires numpy and rasterio: python -m pip install numpy rasterio
Run from any directory; --current defaults to this repository's landcoverchange.
Inputs may be repository roots or landcoverchange directories. Files are paired
by relative path, bands by position (band descriptions must agree when present).
No resampling is performed. Both-NoData pixels are excluded from agreement;
one-sided NoData counts as a difference. Class transitions use jointly valid
pixels only, with the previous reconstruction as the source class.
"""

import argparse
import csv
from collections import Counter
from contextlib import ExitStack
import importlib.util
import os
from pathlib import Path
import sys

# Windows can inherit incompatible PostGIS projection data. Prefer Rasterio's
# bundled data when available, for this Python process only.
if sys.platform == "win32":
    spec = importlib.util.find_spec("rasterio")
    if spec is not None and spec.origin:
        package = Path(spec.origin).parent
        if (package / "proj_data" / "proj.db").is_file():
            os.environ["PROJ_DATA"] = str(package / "proj_data")
            os.environ["PROJ_LIB"] = str(package / "proj_data")
        if (package / "gdal_data").is_dir():
            os.environ["GDAL_DATA"] = str(package / "gdal_data")

import numpy as np
import rasterio
from rasterio.windows import Window


SUMMARY_FIELDS = [
    "file", "band", "age_or_band", "status", "message", "total_pixels",
    "both_nodata", "previous_only_nodata", "current_only_nodata",
    "jointly_valid", "equal_pixels", "class_changed_pixels",
    "different_pixels", "compared_pixels", "agreement_percent",
]
TRANSITION_FIELDS = ["file", "band", "age_or_band", "previous_class", "current_class", "pixels"]


def reconstruction_root(path):
    path = path.expanduser().resolve()
    if not path.is_dir():
        raise ValueError(f"Reconstruction directory does not exist: {path}")
    return path / "landcoverchange" if (path / "landcoverchange").is_dir() else path


def discover(root, pattern):
    return {p.relative_to(root).as_posix(): p for p in root.glob(pattern) if p.is_file()}


def check_alignment(current, previous):
    for attribute in ("crs", "transform", "width", "height", "count"):
        if getattr(current, attribute) != getattr(previous, attribute):
            raise ValueError(f"Different {attribute}: current={getattr(current, attribute)}, "
                             f"previous={getattr(previous, attribute)}")
    for band, (now, old) in enumerate(zip(current.descriptions, previous.descriptions), 1):
        if now and old and now != old:
            raise ValueError(f"Different age/band descriptions at band {band}: {now!r} vs {old!r}")


def compare_pair(current_path, previous_path, relative, summary, transitions, diff_path=None):
    with ExitStack() as stack:
        current = stack.enter_context(rasterio.open(current_path))
        previous = stack.enter_context(rasterio.open(previous_path))
        check_alignment(current, previous)
        difference = None
        if diff_path is not None:
            diff_path.parent.mkdir(parents=True, exist_ok=True)
            # A fresh profile avoids copying data-type-specific TIFF options.
            difference = stack.enter_context(rasterio.open(
                diff_path, "w", driver="GTiff", width=current.width,
                height=current.height, count=current.count, crs=current.crs,
                transform=current.transform, dtype="uint8", nodata=255,
                compress="deflate", tiled=True,
            ))
            difference.update_tags(
                codes="0=equal, 1=class changed, 2=previous NoData, 3=current NoData, 255=both NoData"
            )

        counts_by_band = [Counter() for _ in current.indexes]
        pairs_by_band = [Counter() for _ in current.indexes]
        labels = [now or old or str(band) for band, (now, old) in enumerate(
            zip(current.descriptions, previous.descriptions), 1)]
        for band, label in enumerate(labels, 1):
            if difference is not None:
                difference.set_band_description(band, label)
        # Read each spatial window across all bands together. Terra's stacks
        # interleave bands, so reading one band at a time repeatedly decompresses
        # the same TIFF data. Window size also caps memory for long time series.
        tile_size = min(512, max(1, int((8_000_000 / current.count) ** 0.5)))
        for row in range(0, current.height, tile_size):
            for col in range(0, current.width, tile_size):
                window = Window(col, row, min(tile_size, current.width - col),
                                min(tile_size, current.height - row))
                now_stack = current.read(window=window, masked=True)
                old_stack = previous.read(window=window, masked=True)
                for band in current.indexes:
                    now = now_stack[band - 1]
                    old = old_stack[band - 1]
                    counts = counts_by_band[band - 1]
                    class_pairs = pairs_by_band[band - 1]
                    now_missing = np.ma.getmaskarray(now) | ~np.isfinite(now.data)
                    old_missing = np.ma.getmaskarray(old) | ~np.isfinite(old.data)
                    valid = ~(now_missing | old_missing)
                    both_missing = now_missing & old_missing
                    previous_only = old_missing & ~now_missing
                    current_only = now_missing & ~old_missing
                    changed = valid & (now.data != old.data)
                    counts.update({
                        "both_nodata": int(both_missing.sum()),
                        "previous_only_nodata": int(previous_only.sum()),
                        "current_only_nodata": int(current_only.sum()),
                        "jointly_valid": int(valid.sum()),
                        "class_changed_pixels": int(changed.sum()),
                    })
                    if valid.any() and not changed.any():
                        classes, frequencies = np.unique(now.data[valid], return_counts=True)
                        for value, frequency in zip(classes.tolist(), frequencies):
                            class_pairs[value, value] += int(frequency)
                    elif valid.any():
                        pairs, frequencies = np.unique(
                            np.column_stack((old.data[valid], now.data[valid])),
                            axis=0, return_counts=True,
                        )
                        for pair, frequency in zip(pairs, frequencies):
                            class_pairs[tuple(pair.tolist())] += int(frequency)
                    if difference is not None:
                        codes = np.zeros(now.shape, dtype=np.uint8)
                        codes[changed] = 1
                        codes[previous_only] = 2
                        codes[current_only] = 3
                        codes[both_missing] = 255
                        difference.write(codes, band, window=window)

        any_difference = False
        for band, label in enumerate(labels, 1):
            counts = counts_by_band[band - 1]
            class_pairs = pairs_by_band[band - 1]
            equal = counts["jointly_valid"] - counts["class_changed_pixels"]
            different = (counts["class_changed_pixels"] + counts["previous_only_nodata"]
                         + counts["current_only_nodata"])
            compared = equal + different
            any_difference |= different > 0
            summary.writerow(dict(
                file=relative, band=band, age_or_band=label, status="compared", message="",
                total_pixels=current.width * current.height, **counts,
                equal_pixels=equal, different_pixels=different, compared_pixels=compared,
                agreement_percent=100 * equal / compared if compared else "",
            ))
            for (old_class, now_class), count in sorted(class_pairs.items()):
                transitions.writerow(dict(
                    file=relative, band=band, age_or_band=label,
                    previous_class=old_class, current_class=now_class, pixels=count,
                ))
        return any_difference


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("previous", type=Path, help="Other reconstruction directory")
    parser.add_argument("--current", type=Path,
                        default=Path(__file__).resolve().parents[1] / "landcoverchange")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "landcover_comparison",
                        help="Report directory (existing reports are overwritten)")
    parser.add_argument("--pattern", default="**/LcovOut_new.tif",
                        help="Relative glob in both roots (default: **/LcovOut_new.tif)")
    parser.add_argument("--write-diffs", action="store_true", help="Also save pixel difference TIFFs")
    parser.add_argument("--fail-on-difference", action="store_true", help="Exit 1 when pixels differ")
    args = parser.parse_args(argv)
    try:
        current_root = reconstruction_root(args.current)
        previous_root = reconstruction_root(args.previous)
        current_files = discover(current_root, args.pattern)
        previous_files = discover(previous_root, args.pattern)
        if not current_files or not previous_files:
            raise ValueError(f"No files match {args.pattern!r} in "
                             f"{current_root if not current_files else previous_root}")
    except ValueError as error:
        parser.error(str(error))

    output = args.output.expanduser().resolve()
    # Reports must never overwrite an input, even with a custom glob/output path.
    targets = {output / "summary.csv", output / "transitions.csv"}
    if args.write_diffs:
        targets.update(output / "differences" / name for name in current_files.keys() & previous_files.keys())
    inputs = {p.resolve() for p in (*current_files.values(), *previous_files.values())}
    if any(p.resolve() in inputs for p in targets):
        parser.error("Output paths overlap input rasters; choose another --output directory")
    output.mkdir(parents=True, exist_ok=True)
    issues = 0
    changed_files = 0
    matched = 0
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as summary_file, \
            (output / "transitions.csv").open("w", newline="", encoding="utf-8") as transition_file:
        summary = csv.DictWriter(summary_file, fieldnames=SUMMARY_FIELDS)
        transitions = csv.DictWriter(transition_file, fieldnames=TRANSITION_FIELDS)
        summary.writeheader()
        transitions.writeheader()
        for relative in sorted(current_files.keys() | previous_files.keys()):
            if relative not in current_files or relative not in previous_files:
                status = "missing_current" if relative not in current_files else "missing_previous"
                summary.writerow(dict(file=relative, status=status))
                print(f"{status}: {relative}", file=sys.stderr)
                issues += 1
                continue
            try:
                changed = compare_pair(
                    current_files[relative], previous_files[relative], relative, summary, transitions,
                    output / "differences" / relative if args.write_diffs else None,
                )
                changed_files += int(changed)
                matched += 1
                print(f"{'Different' if changed else 'Identical'}: {relative}")
            except (ValueError, OSError, rasterio.errors.RasterioError) as error:
                summary.writerow(dict(file=relative, status="error", message=str(error)))
                print(f"Cannot compare {relative}: {error}", file=sys.stderr)
                issues += 1
    print(f"Compared {matched} raster pairs; {changed_files} differ; {issues} missing/error pairs.")
    print(f"Reports: {output}")
    return 2 if issues else (1 if args.fail_on_difference and changed_files else 0)


if __name__ == "__main__":
    sys.exit(main())
