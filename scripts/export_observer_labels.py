"""Export PixelFlow class votes by age. Requires NumPy and Rasterio.

From the repository root: python scripts/export_observer_labels.py
Use --ages 3500 -1000 to export just two timeslices.
Site votes are clipped to the fixed output grid.
"""

import argparse
import csv
import importlib.util
import os
import re
import sys
import time
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory

# Select Rasterio's bundled data before importing it on Windows, where an
# inherited PostGIS setting can point to an incompatible projection database.
if sys.platform == "win32":
    spec = importlib.util.find_spec("rasterio")
    if spec is not None and spec.origin:
        package = Path(spec.origin).parent
        proj_data = package / "proj_data"
        if (proj_data / "proj.db").is_file():
            os.environ["PROJ_DATA"] = os.environ["PROJ_LIB"] = str(proj_data)
        if (package / "gdal_data").is_dir():
            os.environ["GDAL_DATA"] = str(package / "gdal_data")

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin
from rasterio.windows import Window


# Fixed band order from tests/observer_labels/band_classes_new.csv.
CLASS_CODES = np.array([50, 65, 75, 85, 90, 105, 125, 130, 140, 155, 165,
                        180, 205, 210, 220])
MASK_BAND = len(CLASS_CODES) + 1
TILE_SIZE = 512


def target_grid(reference):
    if reference is None:
        return CRS.from_epsg(4326), from_origin(-172, 72, 1 / 360, 1 / 360), 18000, 8280
    with rasterio.open(reference) as dataset:
        transform = dataset.transform
        if (dataset.crs is None or transform.b != 0 or transform.d != 0
                or transform.a <= 0 or transform.e >= 0):
            raise ValueError("The grid reference must have a CRS and be north-up.")
        return dataset.crs, transform, dataset.width, dataset.height


def site_offset(dataset, grid):
    crs, transform, _, _ = grid
    # Check both corners in pixel coordinates to catch origin and scale drift.
    relative = ~transform * dataset.transform
    col, row = (round(value) for value in relative * (0, 0))
    for x, y in ((0, 0), (dataset.width, 0), (0, dataset.height)):
        if not np.allclose(relative * (x, y), (col + x, row + y), rtol=0, atol=1e-6):
            raise ValueError(f"Site pixels do not align to the output grid: {dataset.name}")
    if dataset.crs != crs:
        raise ValueError(f"Site CRS differs from the output grid: {dataset.name}")
    return col, row


def index_timeslices(source, grid):
    status_file = source / "ice_flow_status.csv"
    with status_file.open(newline="", encoding="utf-8-sig") as stream:
        status = csv.DictReader(stream)
        if not {"ID", "status"}.issubset(status.fieldnames or []):
            raise ValueError(f"Expected ID and status columns in {status_file}")
        completed = {int(row["ID"]) for row in status if row["status"] in {"ok", "cached"}}

    site_paths = {}
    for path in sorted(source.glob("PSA_*/psaChange/LcovOut_new.tif")):
        match = re.fullmatch(r"PSA_(\d+)(?:__.+)?", path.parent.parent.name)
        if match is None or int(match[1]) not in completed:
            continue
        site = int(match[1])
        if site in site_paths:
            raise ValueError(f"Multiple completed TIFFs for PSA {site}: {site_paths[site]}, {path}")
        site_paths[site] = path
    missing = sorted(completed - site_paths.keys())
    if missing:
        raise ValueError(f"Missing LcovOut_new.tif for completed sites: {missing}")

    timeslices = defaultdict(list)
    for path in site_paths.values():
        with rasterio.open(path) as dataset:
            col, row = site_offset(dataset, grid)
            footprint = Window(col, row, dataset.width, dataset.height)
            ages = set()
            for band, name in enumerate(dataset.descriptions, 1):
                match = re.fullmatch(r"LCov_(-?\d+)", name or "")
                if match is None:
                    raise ValueError(f"Expected an LCov_<age> description on band {band}: {path}")
                age = int(match[1])
                if age in ages:
                    raise ValueError(f"Duplicate LCov_{age} band: {path}")
                ages.add(age)
                timeslices[age].append((path, band, footprint))
    if not timeslices:
        raise ValueError(f"No completed PixelFlow timeslices found in {source}")
    return timeslices


def add_votes(votes, tile, dataset, band, footprint):
    left = max(tile.col_off, footprint.col_off)
    top = max(tile.row_off, footprint.row_off)
    right = min(tile.col_off + tile.width, footprint.col_off + footprint.width)
    bottom = min(tile.row_off + tile.height, footprint.row_off + footprint.height)
    source_window = Window(left - footprint.col_off, top - footprint.row_off,
                           right - left, bottom - top)
    labels = dataset.read(band, window=source_window, masked=True)
    rows, cols = np.nonzero(~np.ma.getmaskarray(labels))
    values = labels.data[rows, cols]
    indices = np.searchsorted(CLASS_CODES, values)
    recognized = CLASS_CODES[indices.clip(max=len(CLASS_CODES) - 1)] == values
    if not recognized.all():
        unknown = np.unique(values[~recognized]).tolist()
        raise ValueError(f"Unmapped classes {unknown} in {dataset.name}, band {band}")
    votes[indices, rows + top - tile.row_off, cols + left - tile.col_off] += 1


def publish(staging, destination):
    # Windows readers such as antivirus scanners can briefly hold a TIFF open.
    for delay in (0.1, 0.25, 0.5, 1):
        try:
            staging.rename(destination)
            return
        except PermissionError as error:
            if getattr(error, "winerror", None) not in {5, 32, 33} or destination.exists():
                raise
            time.sleep(delay)
    staging.rename(destination)


def export_timeslice(destination, inputs, grid):
    crs, transform, width, height = grid
    # Publish the directory only after both files have been written and closed.
    with TemporaryDirectory(dir=destination.parent, prefix=f".{destination.name}_") as temp:
        staging = Path(temp) / destination.name
        staging.mkdir()
        with ExitStack() as stack:
            tile_sites = defaultdict(list)
            for path, band, footprint in inputs:
                # Clip the tile range, keeping the original source origin for reads.
                left = max(0, footprint.col_off)
                top = max(0, footprint.row_off)
                right = min(width, footprint.col_off + footprint.width)
                bottom = min(height, footprint.row_off + footprint.height)
                if left >= right or top >= bottom:
                    continue
                dataset = stack.enter_context(rasterio.open(path))
                for row in range(top // TILE_SIZE, (bottom - 1) // TILE_SIZE + 1):
                    for col in range(left // TILE_SIZE, (right - 1) // TILE_SIZE + 1):
                        tile_sites[row, col].append((dataset, band, footprint))
            output = stack.enter_context(rasterio.open(
                staging / "vote_counts_new.tif", "w", driver="GTiff",
                crs=crs, transform=transform, width=width, height=height,
                count=MASK_BAND, dtype="uint8", compress="deflate", tiled=True,
                blockxsize=TILE_SIZE, blockysize=TILE_SIZE, BIGTIFF="IF_SAFER",
            ))
            # GeoTIFF fills unwritten tiles with zero when there is no NoData value.
            for (row, col), sites in sorted(tile_sites.items()):
                tile = output.block_window(1, row, col)
                max_vote = np.iinfo(np.uint8).max
                dtype = np.uint8 if len(sites) <= max_vote else np.uint32
                votes = np.zeros((MASK_BAND, tile.height, tile.width), dtype=dtype)
                for dataset, band, footprint in sites:
                    add_votes(votes, tile, dataset, band, footprint)
                if dtype == np.uint32 and votes[:-1].max() > max_vote:
                    raise ValueError(f"A class vote count exceeds 255 in {destination.name}, {tile}")
                votes[-1] = np.any(votes[:-1], axis=0)
                output.write(votes.astype(np.uint8, copy=False), window=tile)

        with (staging / "band_classes_new.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("band_index", "class_code"))
            writer.writerows((band, f"{code}.0") for band, code in enumerate(CLASS_CODES, 1))
            writer.writerow((MASK_BAND, "mask"))
        publish(staging, destination)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path,
                        default=Path(__file__).resolve().parents[1] / "landcoverchange",
                        help="PixelFlow results directory containing ice_flow_status.csv")
    parser.add_argument("--output", type=Path, help="Output root (default: SOURCE/observer_labels)")
    parser.add_argument("--ages", type=int, nargs="+", help="Exact ages to export (default: all)")
    parser.add_argument("--grid-tif", type=Path,
                        help="Aligned grid template (default: the example's Alaska grid)")
    args = parser.parse_args(argv)

    try:
        grid = target_grid(args.grid_tif)
        timeslices = index_timeslices(args.source, grid)
        ages = sorted(set(args.ages) if args.ages else timeslices)
        missing = sorted(set(ages) - timeslices.keys())
        if missing:
            raise ValueError(f"No completed sites have these ages: {missing}")
        root = args.output or args.source / "observer_labels"
        destinations = {age: root / f"LCov_{age}" for age in ages}
        existing = [str(path) for path in destinations.values() if path.exists()]
        if existing:
            raise ValueError(f"Output directories already exist: {', '.join(existing)}")
        root.mkdir(parents=True, exist_ok=True)
        for age, destination in destinations.items():
            export_timeslice(destination, timeslices[age], grid)
            print(f"{destination} ({len(timeslices[age])} sites)", flush=True)
    except (OSError, ValueError, csv.Error, rasterio.errors.RasterioError) as error:
        parser.exit(1, f"Export failed: {error}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
