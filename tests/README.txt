Pixel comparison of reconstruction versions
==========================================
Install dependencies in your Python environment:
  python -m pip install numpy rasterio

From the repository root:
  python tests/compare_landcover.py "C:/path/to/previous/landcoverchange"

The current reconstruction defaults to this repository's landcoverchange folder.
Both inputs also accept a repository root containing a landcoverchange folder.
For example:
  python tests/compare_landcover.py "D:/previous" --current "D:/current" --output "tests/comparison_v1_v2" --write-diffs

Final LcovOut_new.tif stacks are matched by relative path (site directory).
Every band and pixel is compared exactly. Dimensions, CRS, transform and band
counts must agree. Band descriptions, when present in both files, must agree
in the same order. Bands without descriptions are compared by position; check
their age order before using such files. The script never resamples rasters.
To select other final filenames or individual age TIFFs in both versions, use
--pattern "**/LCov_*.tif" (or another relative glob).

Outputs (overwritten on rerun):
  summary.csv: per-file/per-band counts, agreement and missing/error records.
  transitions.csv: previous_class -> current_class pixel counts, including
    unchanged classes, for pixels valid in both reconstructions.
  differences/: optional georeferenced TIFFs with the same bands as the input:
    0 = equal, 1 = changed class, 2 = previous NoData, 3 = current NoData,
    255 = both NoData (output NoData value).

Both-NoData pixels are excluded from agreement. Pixels valid in only one
version count as differences. NaN/infinite values are treated as NoData.
Agreement is 100 * equal / (equal + class changes + one-sided NoData).
Agreement is empty if all pixels in a band are NoData in both versions.

Exit codes: 0 = comparison completed; 2 = missing files or comparison errors.
With --fail-on-difference, completed comparisons with differences exit 1.
Missing/error pairs still exit 2. Differences are expected for changed models;
the flag is optional for regression checks. Read summary.csv for incomplete
comparisons, including a possible error after earlier bands were written.
On Windows, Rasterio's bundled projection/GDAL data are selected when present
to avoid conflicts with inherited PostGIS settings, for this process only.

Run synthetic integration checks:
  python -m unittest discover -s tests -p test_compare_landcover.py
