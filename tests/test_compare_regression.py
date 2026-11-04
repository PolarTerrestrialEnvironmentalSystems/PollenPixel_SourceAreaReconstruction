"""Tests for saved-score comparison, including R serialization/extraction."""

import contextlib
import csv
import io
from pathlib import Path
import subprocess
import tempfile
import unittest

from compare_regression import compare_rows, find_rscript, main


def metric(version, value="0.2", **overrides):
    row = dict(version=version, file="site/summaryResults/veg_lc_model_new.rda",
               model="calibMap_1", scope="overall", **{"class": ""}, metric="mean_MAE",
               value=value, direction="lower", status="ok", method="ElasticNet_CLR",
               class_set="65;75", samples="10", folds="2", pseudocount="1e-06",
               calibration_digest="data", fold_digest="folds")
    row.update(overrides)
    return row


class ScoreTests(unittest.TestCase):
    def test_directions_and_percent_change(self):
        rows, errors = compare_rows([metric("previous"), metric("current", "0.1")])
        self.assertFalse(errors)
        self.assertAlmostEqual(rows[0]["delta"], -0.1)
        self.assertAlmostEqual(rows[0]["percent_change"], -50)
        self.assertEqual(rows[0]["assessment"], "improved")
        rows, _ = compare_rows([metric("previous", "-0.5", metric="mean_R2", direction="higher"),
                                metric("current", "-0.2", metric="mean_R2", direction="higher")])
        self.assertEqual(rows[0]["assessment"], "improved")
        self.assertAlmostEqual(rows[0]["percent_change"], 60)

    def test_changed_or_unknown_context_prevents_claim(self):
        for overrides in ({"class_set": "65;75;140"}, {"fold_digest": "new"},
                          {"calibration_digest": ""}):
            rows, _ = compare_rows([metric("previous"), metric("current", "0.1", **overrides)])
            self.assertEqual(rows[0]["assessment"], "context_requires_review")

    def test_unavailable_zero_and_missing_metrics(self):
        for value in ("", "NaN", "Inf"):
            rows, _ = compare_rows([metric("previous", value), metric("current")])
            self.assertEqual(rows[0]["status"], "unavailable_value")
            self.assertNotIn("delta", rows[0])
        rows, _ = compare_rows([metric("previous", "0"), metric("current")])
        self.assertEqual(rows[0]["percent_change"], "")
        rows, _ = compare_rows([metric("previous", metric="legacy_mean_MAE"), metric("current")])
        self.assertEqual({r["status"] for r in rows}, {"missing_previous_metric", "missing_current_metric"})

    def test_duplicate_metric_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            compare_rows([metric("current"), metric("current")])


class ExtractionTests(unittest.TestCase):
    def test_serialized_models_and_missing_file(self):
        try:
            rscript = find_rscript(None)
        except ValueError as error:
            self.skipTest(str(error))
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            root = Path(directory).resolve()
            fixture = root / "fixture.R"
            fixture.write_text('''
root <- commandArgs(trailingOnly=TRUE)[1]
for (version in c("current", "previous")) {
  value <- if (version == "current") 0.1 else 0.2
  m <- list(overall_metrics=data.frame(method="ElasticNet_CLR",mean_MAE=value),
            per_class_metrics=data.frame(class=c("65","75"),MAE=c(value,value)),
            samplewise_metrics=data.frame(fold=c(1L,2L)),pseudocount=1e-6,
            model_object=list(calibration_data=data.frame(Dataset_ID=c(1L,2L),calibrationMap=1L)))
  legacy <- list(mean_mae=5,mean_mae_per_class=c("65"=4,"75"=6))
  veg_lc_model <- list(calibMap_1=m,calibMap_2=legacy)
  folder <- file.path(root,version,"site","summaryResults")
  dir.create(folder,recursive=TRUE)
  save(veg_lc_model,file=file.path(folder,"veg_lc_model_new.rda"))
}
''', encoding="utf-8")
            subprocess.run([rscript, "--vanilla", str(fixture), str(root)],
                           check=True, capture_output=True)
            arguments = [str(root / "previous"), "--current", str(root / "current"),
                         "--output", str(root / "report"), "--rscript", rscript]
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                code = main(arguments)
            self.assertEqual(code, 0)
            with (root / "report/metrics.csv").open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 6)
            modern = [r for r in rows if not r["metric"].startswith("legacy_")]
            self.assertTrue(all(r["assessment"] == "improved" for r in modern))
            (root / "current/site/summaryResults/veg_lc_model_new.rda").rename(
                root / "current/site/summaryResults/other.rda")
            extra = root / "current/extra/summaryResults"
            extra.mkdir(parents=True)
            (extra / "veg_lc_model_new.rda").write_bytes(b"invalid R file")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                code = main(arguments)
            self.assertEqual(code, 2)
            with (root / "report/issues.csv").open(newline="", encoding="utf-8") as stream:
                issues = list(csv.DictReader(stream))
            self.assertEqual(len(issues), 3)  # invalid file plus two missing counterparts


if __name__ == "__main__":
    unittest.main()
