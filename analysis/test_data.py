#!/usr/bin/env python3
"""Independent scalar checks; optional reviewed upstream CLI cross-check.

python test_data.py --package package --casland Casland
"""
import argparse
import csv
import io
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--package", type=Path, default=Path(__file__).parent / "package")
    ap.add_argument("--casland", type=Path)
    args = ap.parse_args()
    data = json.loads((args.package / "data.json").read_text())
    guides = data["guides"]
    assert len(guides) == 91
    n, missing, calc_review = 0, 0, 0
    for guide in guides:
        ts, cm, ko = guide["times"], guide["controls"]["CM"], guide["controls"]["KO"]
        assert len(ts) == (15 if guide["id"] == "S32" else 16)
        assert len(guide["targets"]) == 61
        zero = [i for i, t in enumerate(ts) if t == 0]
        def late_fraction(counts, control, keep):
            baseline = sum(counts[i] for i in keep)
            control0 = sum(control[i] for i in keep)
            if baseline == 0 or control0 == 0:
                return None
            return sum(1 - counts[i] / (baseline / control0 * control[i] + .1) for i in [-2, -1]) / 2
        for target in guide["targets"]:
            counts = target["counts"]
            n += len(ts)
            pooled = [a + b for a, b in zip(cm, ko)]
            baseline = sum(counts[i] for i in zero)
            assert baseline == target["baseline_reads"]
            value = late_fraction(counts, pooled, zero)
            if value is None:
                missing += 1
                assert target["late"] is None
                assert all(x is None for x in target["f_bound"])
            else:
                assert abs(value - target["late"]) <= .00000051
                shift = abs(late_fraction(counts, cm, zero) - late_fraction(counts, ko, zero))
                assert abs(shift - target["control_shift"]) <= .00000051
                alts = [late_fraction(counts, pooled, [j for j in zero if j != i]) for i in zero]
                if any(x is None for x in alts):
                    assert target["baseline_shift"] is None
                else:
                    assert abs(max(alts) - min(alts) - target["baseline_shift"]) <= .00000051
                gap = abs(target["f_bound"][-2] - target["f_bound"][-1])
                assert abs(gap - target["late_pair_gap"]) <= .00000101
            if baseline < 100 or any(target[k] is None or target[k] > .15 for k in ["control_shift", "baseline_shift", "late_pair_gap"]):
                calc_review += 1
    assert n == 88755
    assert calc_review == 1514
    csv_path = args.package / "drilldown_counts.csv"
    opener = csv_path.open if csv_path.exists() else lambda: gzip.open(str(csv_path) + ".gz", "rt")
    with opener() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == n
    lookup = {(g["id"], t["id"]): (g, t) for g in guides for t in g["targets"]}
    for row in rows:
        g, t = lookup[row["sublibrary"], row["element"]]
        i = g["timepoint_ids"].index(int(row["timepoint_id"]))
        assert float(row["time_min"]) == g["times"][i]
        assert int(row["count"]) == t["counts"][i]
    report = {"scalar_normalization_checked_targets": 5551, "csv_raw_count_rows_checked": n,
              "zero_baseline_targets": missing, "default_review_count_checked": calc_review}
    if args.casland:
        # Execution is limited to the exact reviewed upstream script and table.
        for name, expected in {"list_ddG_perturbations.py": "8ef9769123bcddf8a8da999c154966612313342a44a363b9aeb5698b14320e62",
                               "S17_F7C_ddG_perturbation.tsv": "01d93415551daa556a81166e77df73eb6e7eaca3f4b9626fa9aff284406f2a97"}.items():
            assert hashlib.sha256((args.casland / name).read_bytes()).hexdigest() == expected
        guide = next(g for g in guides if g["id"] == "S60")
        result = subprocess.run([sys.executable, str(args.casland / "list_ddG_perturbations.py"),
                                 "--ddG_table", str(args.casland / "S17_F7C_ddG_perturbation.tsv"),
                                 guide["sequence"], "binding"], check=True, capture_output=True, text=True)
        model = {r["target_name"]: r for r in csv.DictReader(io.StringIO(result.stdout), delimiter="\t")}
        tested = 0
        for target in guide["targets"]:
            if target["mutation"] == "WT":
                continue
            record = model[target["model_name"]]
            assert record["target_sequence"] == target["sequence"], target
            assert abs(float(record["ddG_perturbation"]) - target["model_ddg"]) <= .00000051
            tested += 1
        report["upstream_cli_L1_single_variants_checked"] = tested
        report["upstream_commit"] = "9193bf679c2bab346d8c0f30a1f92e52c29c9bed"
        (args.package / "casland_L1_binding_cli.tsv").write_text(result.stdout)
    (args.package / "validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
