"""Synthetic integration checks: python -m unittest discover -s tests -p test_compare_landcover.py"""

import contextlib
import csv
import io
from pathlib import Path
import tempfile
import unittest

from compare_landcover import main

import numpy as np
import rasterio
from rasterio.transform import from_origin


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.current = self.root / "current"
        self.previous = self.root / "previous"
        self.output = self.root / "report"
        self.relative = Path("PSA_1/psaChange/LcovOut_new.tif")

    def write(self, root, values, *, transform=None, descriptions=None):
        path = root / self.relative
        path.parent.mkdir(parents=True, exist_ok=True)
        data = np.asarray(values, dtype="int16")
        with rasterio.open(
            path, "w", driver="GTiff", width=data.shape[2], height=data.shape[1],
            count=data.shape[0], dtype=data.dtype, nodata=-9999, crs="EPSG:4326",
            transform=transform or from_origin(0, 10, 1, 1),
        ) as dataset:
            dataset.write(data)
            dataset.descriptions = descriptions or tuple(f"LCov_{i * 500}" for i in range(len(data)))

    def run_comparison(self, *options):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = main([str(self.previous), "--current", str(self.current),
                         "--output", str(self.output), *options])
        with (self.output / "summary.csv").open(newline="", encoding="utf-8") as report:
            return code, list(csv.DictReader(report))

    def test_changes_masks_transitions_and_difference_raster(self):
        self.write(self.previous, [[[10, 20, -9999], [30, -9999, 210]]])
        self.write(self.current, [[[10, 40, 50], [-9999, -9999, 210]]])
        code, rows = self.run_comparison("--write-diffs", "--fail-on-difference")
        self.assertEqual(code, 1)
        row = rows[0]
        for field, value in {"both_nodata": 1, "previous_only_nodata": 1,
                             "current_only_nodata": 1, "class_changed_pixels": 1,
                             "equal_pixels": 2, "different_pixels": 3,
                             "compared_pixels": 5, "agreement_percent": 40}.items():
            self.assertEqual(float(row[field]), value, field)
        with rasterio.open(self.output / "differences" / self.relative) as dataset:
            np.testing.assert_array_equal(dataset.read(1), [[0, 1, 2], [3, 255, 0]])
            self.assertEqual(dataset.descriptions, ("LCov_0",))
        with (self.output / "transitions.csv").open(newline="", encoding="utf-8") as report:
            pairs = {(float(r["previous_class"]), float(r["current_class"])): int(r["pixels"])
                     for r in csv.DictReader(report)}
        self.assertEqual(pairs, {(10, 10): 1, (20, 40): 1, (210, 210): 1})

    def test_identical_multiband_and_all_nodata(self):
        values = [[[10, 220]], [[-9999, -9999]]]
        self.write(self.current, values)
        self.write(self.previous, values)
        code, rows = self.run_comparison("--fail-on-difference")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 2)
        self.assertEqual(float(rows[0]["agreement_percent"]), 100)
        self.assertEqual(rows[1]["agreement_percent"], "")

    def test_shifted_grid_is_rejected(self):
        self.write(self.current, [[[10]]])
        self.write(self.previous, [[[10]]], transform=from_origin(1, 10, 1, 1))
        code, rows = self.run_comparison()
        self.assertEqual(code, 2)
        self.assertEqual(rows[0]["status"], "error")
        self.assertIn("transform", rows[0]["message"])

    def test_band_order_is_rejected(self):
        self.write(self.current, [[[10]], [[20]]])
        self.write(self.previous, [[[10]], [[20]]], descriptions=("LCov_500", "LCov_0"))
        code, rows = self.run_comparison()
        self.assertEqual(code, 2)
        self.assertIn("descriptions", rows[0]["message"])

    def test_missing_files_in_either_version_are_reported(self):
        self.write(self.current, [[[10]]])
        self.relative = Path("PSA_2/psaChange/LcovOut_new.tif")
        self.write(self.previous, [[[10]]])
        code, rows = self.run_comparison()
        self.assertEqual(code, 2)
        self.assertEqual({r["status"] for r in rows}, {"missing_current", "missing_previous"})


if __name__ == "__main__":
    unittest.main()
