"""Write a single-band majority_vote.tif into each observer-label time slice.

Requires NumPy and Rasterio. From the repository root:
  python scripts/export_majority_vote.py

The highest vote count wins. For each tie, a class is chosen uniformly
at random from the classes sharing the highest vote count.
Pixels outside the coverage mask or without votes are NoData (0).
"""

import argparse
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

# Reuse the observer exporter's tile size and Windows GDAL/PROJ setup.
from export_observer_labels import TILE_SIZE, np, rasterio


def read_band_mapping(path, band_count):
    bands = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not {"band_index", "class_code"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Expected band_index and class_code columns in {path}")
        for row in reader:
            band = int(row["band_index"] or "")
            label = (row["class_code"] or "").strip()
            if band in bands:
                raise ValueError(f"Duplicate band index {band} in {path}")
            if label == "mask":
                bands[band] = None
            else:
                code = float(label)
                if not code.is_integer() or not 1 <= code <= 255:
                    raise ValueError(f"Expected an integer class code from 1 to 255 in {path}")
                bands[band] = int(code)

    mask_bands = [band for band, code in bands.items() if code is None]
    classes = sorted((code, band) for band, code in bands.items() if code is not None)
    if set(bands) != set(range(1, band_count + 1)) or len(mask_bands) != 1 or not classes:
        raise ValueError(f"Expected one mapping per raster band, including one mask, in {path}")
    codes, class_bands = zip(*classes)
    if len(set(codes)) != len(codes):
        raise ValueError(f"Duplicate class codes in {path}")
    return np.array(codes, dtype=np.uint8), [*class_bands, mask_bands[0]]


def export_majority_vote(directory):
    destination = directory / "majority_vote.tif"
    rng = np.random.default_rng()
    with rasterio.open(directory / "vote_counts_new.tif") as source:
        codes, bands = read_band_mapping(directory / "band_classes_new.csv", source.count)
        profile = source.profile
        profile.update(driver="GTiff", count=1, dtype="uint8", nodata=0,
                       compress="deflate", tiled=True,
                       blockxsize=TILE_SIZE, blockysize=TILE_SIZE)

        # Close the completed TIFF before replacing the destination on reruns.
        with TemporaryDirectory(dir=directory, prefix=".majority_vote_") as temp:
            staging = Path(temp) / destination.name
            with rasterio.open(staging, "w", **profile) as output:
                output.set_band_description(1, "land_cover_class")
                for _, window in output.block_windows(1):
                    data = source.read(bands, window=window)
                    votes, coverage = data[:-1], data[-1]
                    winning_votes = votes.max(axis=0)
                    valid = (coverage != 0) & (winning_votes != 0)
                    top_classes = votes == winning_votes
                    winners = top_classes.argmax(axis=0)
                    top_counts = top_classes.sum(axis=0)
                    tied = valid & (top_counts > 1)
                    if tied.any():
                        # Sample a rank among each pixel's tied classes.
                        ranks = rng.integers(top_counts[tied])
                        cumulative = top_classes[:, tied].cumsum(axis=0)
                        winners[tied] = (cumulative > ranks).argmax(axis=0)
                    labels = codes[winners]
                    labels[~valid] = 0
                    output.write(labels, 1, window=window)
            staging.replace(destination)
    return destination


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path,
                        default=Path(__file__).resolve().parents[1] / "landcoverchange" / "observer_labels",
                        help="Root containing the time-slice folders (default: landcoverchange/observer_labels)")
    args = parser.parse_args(argv)

    try:
        directories = sorted(path for path in args.source.iterdir() if path.is_dir() and (
            (path / "vote_counts_new.tif").is_file() or (path / "band_classes_new.csv").is_file()
        ))
        if not directories:
            raise ValueError(f"No observer-label time slices found in {args.source}")
        for directory in directories:
            print(export_majority_vote(directory), flush=True)
    except (OSError, ValueError, csv.Error, rasterio.errors.RasterioError) as error:
        parser.exit(1, f"Majority vote export failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
