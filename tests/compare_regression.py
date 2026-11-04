"""Compare saved cross-validation metrics; requires Python and R, no Python packages.

Inputs follow compare_landcover.py: repository, landcoverchange, or matching
site directories. Models are paired by relative file path and calibration-map
name. No fitting is performed. Deltas are current minus previous. Modern
MAE/RMSE values are proportions (0-1), not percentage points. Results describe
saved CV scores, not statistical significance or an independent test set.
"""

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


KEY = ("file", "model", "scope", "class", "metric")
CONTEXT = ("method", "class_set", "samples", "folds", "pseudocount",
           "calibration_digest", "fold_digest")
FIELDS = [*KEY, "status", "previous", "current", "delta", "percent_change",
          "better_direction", "assessment", "context", "context_details"]


def root(path):
    path = path.expanduser().resolve()
    if not path.is_dir():
        raise ValueError(f"Directory does not exist: {path}")
    return path / "landcoverchange" if (path / "landcoverchange").is_dir() else path


def find_rscript(explicit):
    if explicit:
        return explicit
    found = shutil.which("Rscript")
    if found:
        return found
    if os.name == "nt":
        candidates = list((Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "R").glob("R-*/bin/Rscript.exe"))
        def version(path):
            return tuple(int(p) for p in path.parents[1].name.removeprefix("R-").split(".") if p.isdigit())
        if candidates:
            return str(max(candidates, key=version))
    raise ValueError("Rscript was not found. Install R or provide --rscript PATH.")


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def compare_rows(rows):
    indexed = {"previous": {}, "current": {}}
    errors = []
    for row in rows:
        if row["status"] != "ok":
            errors.append(row)
            continue
        key = tuple(row[k] for k in KEY)
        table = indexed[row["version"]]
        if key in table:
            raise ValueError(f"Duplicate metric key: {key}")
        table[key] = row
    result = []
    for key in sorted(indexed["previous"].keys() | indexed["current"].keys()):
        old = indexed["previous"].get(key)
        now = indexed["current"].get(key)
        row = dict(zip(KEY, key))
        row.update(previous=old["value"] if old else "", current=now["value"] if now else "")
        if not old or not now:
            row["status"] = "missing_previous_metric" if not old else "missing_current_metric"
            result.append(row)
            continue
        changed = [field for field in CONTEXT if old.get(field) != now.get(field)]
        unknown = [field for field in CONTEXT if not old.get(field) or not now.get(field)]
        context = "changed" if changed else "unknown" if unknown else "same_recorded_context"
        row.update(context=context, context_details="; ".join(
            [f"{field}: {old.get(field, '')} -> {now.get(field, '')}" for field in changed]
            + [f"{field}: unavailable" for field in unknown if field not in changed]))
        old_value, now_value = number(old["value"]), number(now["value"])
        row["better_direction"] = now["direction"]
        if old_value is None or now_value is None:
            row.update(status="unavailable_value", assessment="unavailable")
        else:
            delta = now_value - old_value
            row.update(status="compared", delta=delta,
                       percent_change=100 * delta / abs(old_value) if old_value else "")
            if delta == 0:
                row["assessment"] = "unchanged"
            elif context != "same_recorded_context":
                row["assessment"] = "context_requires_review"
            else:
                improved = delta > 0 if now["direction"] == "higher" else delta < 0
                row["assessment"] = "improved" if improved else "regressed"
        result.append(row)
    return result, errors


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("previous", type=Path)
    parser.add_argument("--current", type=Path, default=Path(__file__).resolve().parents[1] / "landcoverchange")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "regression_comparison")
    parser.add_argument("--pattern", default="**/veg_lc_model_new.rda", help="Relative glob in both versions")
    parser.add_argument("--rscript", help="Path to Rscript executable (otherwise detected)")
    args = parser.parse_args(argv)
    try:
        roots = {"previous": root(args.previous), "current": root(args.current)}
        executable = find_rscript(args.rscript)
        manifest = []
        files = {}
        for version, directory in roots.items():
            files[version] = {p.relative_to(directory).as_posix(): p for p in directory.glob(args.pattern) if p.is_file()}
            if not files[version]:
                raise ValueError(f"No files match {args.pattern!r} in {directory}")
            manifest.extend(dict(version=version, file=name, path=str(path)) for name, path in sorted(files[version].items()))
        output = args.output.expanduser().resolve()
        inputs = {p.resolve() for version in files.values() for p in version.values()}
        if any((output / name).resolve() in inputs for name in ("metrics.csv", "saved_metrics.csv", "issues.csv")):
            raise ValueError("Output paths overlap input files; choose another --output directory")
        output.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="regression_", dir=output) as temp:
            manifest_path = Path(temp) / "manifest.csv"
            extracted = Path(temp) / "metrics.csv"
            write_csv(manifest_path, ["version", "file", "path"], manifest)
            process = subprocess.run([executable, "--vanilla", str(Path(__file__).with_name("export_regression_metrics.R")),
                                      str(manifest_path), str(extracted)], capture_output=True, text=True, encoding="utf-8", errors="replace")
            if process.returncode:
                raise ValueError(f"R extraction failed:\n{process.stderr}")
            if process.stderr.strip():
                print(process.stderr.strip(), file=sys.stderr)
            with extracted.open(encoding="utf-8-sig", newline="") as stream:
                reader = csv.DictReader(stream)
                raw_fields = reader.fieldnames
                rows = list(reader)
            comparisons, extraction_errors = compare_rows(rows)
            write_csv(output / "saved_metrics.csv", raw_fields, rows)
        issues = [dict(version=r["version"], file=r["file"], message=r["message"]) for r in extraction_errors]
        for version, other in (("previous", "current"), ("current", "previous")):
            issues.extend(dict(version=other, file=name, message="Missing model file")
                          for name in sorted(files[version].keys() - files[other].keys()))
        write_csv(output / "metrics.csv", FIELDS, comparisons)
        write_csv(output / "issues.csv", ["version", "file", "message"], issues)
        compared = sum(row["status"] == "compared" for row in comparisons)
        unmatched = sum(row["status"].startswith("missing_") for row in comparisons)
        print(f"Compared {compared} metrics; {unmatched} unmatched metrics; {len(issues)} file issues.")
        print(f"Reports: {output}")
        return 2 if issues or not compared else 0
    except (ValueError, OSError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    sys.exit(main())
