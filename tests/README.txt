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

Regression performance comparison
=================================
Requires Python (standard library only) and R; no model refitting or extra R
packages are needed to extract saved metrics. Rscript is detected on PATH or
in the standard Windows R installation directory; override with --rscript.

  python tests/compare_regression.py "D:/previous/landcoverchange"
  python tests/compare_regression.py "D:/previous/site" --current "D:/current/site" --output tests/regression_v1_v2

As with the pixel comparison, pass matching directory levels. Model files are
matched by relative path; calibration maps by their saved names. The default
pattern is **/veg_lc_model_new.rda; --pattern changes it for both inputs.

Outputs (overwritten on rerun):
  metrics.csv: overall and per-class cross-validation scores, current minus
    previous delta, percent change, preferred direction and context checks.
  saved_metrics.csv: extracted metrics and context from each version.
  issues.csv: missing files and unsupported/unreadable model files.

Modern metrics include MAE, RMSE, R2 and mean Aitchison distance where saved.
MAE/RMSE are in proportions (0-1), not percentages. Lower errors/distances and
higher R2 are better. Percent change uses abs(previous) as denominator; it is
blank for a zero baseline. NA/infinite values remain unavailable, not zero.

The comparison checks saved method, class set, sample/fold counts, pseudocount,
calibration row identities/order and fold assignments. Missing or changed
context prevents an improved/regressed assessment for changed scores. Matching
context does not prove identical input values or training implementations.
Review context_details, and review unmatched metrics/classes. Class additions
affect overall averages and compositional distances even for unchanged classes.
No significance testing or paired sample comparison is attempted. The reported
scores are the workflow's saved cross-validation estimates, not an independent
test set; model selection protocols may differ between versions.

Legacy mean_MAE/MAE metrics retain legacy_ prefixes and their stored units.
They are not equated with modern composition metrics. Unsupported legacy
formats appear in issues.csv; missing metrics appear in metrics.csv. Available
metrics are still reported when other files fail. Exit 2 means file issues or
no numeric metric pairs; exit 0 means the comparison completed (inspect the
status/context columns for unavailable or unmatched metrics).

Run regression comparison checks:
  python -m unittest discover -s tests -p test_compare_regression.py

Worker configuration checks (no production data modified):
  Rscript tests/workers.R
