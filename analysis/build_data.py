#!/usr/bin/env python3
"""Independent sensitivity audit of deposited filter-binding counts, not the authors' fits.

Run: python build_data.py --source . --output package
Requires numpy, pandas, openpyxl. Does not execute downloaded code or use network.
"""
import argparse
import csv
import hashlib
import json
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

COMPLEMENT = str.maketrans("ATGC", "TACG")
TRANSITION = str.maketrans("ATGC", "GCAT")
OPPOSING = str.maketrans("ATGC", "CGTA")
SINGLE = re.compile(r"^-(?:[1-9]|1[0-9]|20):[ACGT]$")
COMMIT = "9193bf679c2bab346d8c0f30a1f92e52c29c9bed"
SOURCE_MD5 = {"count_data.zip": "f123c9defb2917480f9b4483e44bd43e", "HiTBiP_libraries.xlsx": "4cb9e86dc6da79e9f643753f514a8ad7"}
CASLAND_SHA256 = {"S17_F7C_ddG_perturbation.tsv": "01d93415551daa556a81166e77df73eb6e7eaca3f4b9626fa9aff284406f2a97",
                  "list_ddG_perturbations.py": "8ef9769123bcddf8a8da999c154966612313342a44a363b9aeb5698b14320e62",
                  "LICENSE": "591b424ba29232a92a4da8234dac3a7404de75a55d28da640ae776e4149d0621"}


def floats(x):
    return [round(float(v), 6) if np.isfinite(v) else None for v in np.asarray(x).ravel()]


def number(x):
    return round(float(x), 6) if np.isfinite(x) else None


def normalized(counts, controls, zero):
    """Pool metadata-listed zero-time points as exposure-normalized count rate.

    Expected = sum(target zero counts) / sum(control zero counts) * control_t.
    Fraction depleted = 1 - observed / (expected + 0.1), per paper plotting equation.
    A zero target baseline is unidentifiable (NaN), not an apparent negative binder.
    """
    target0 = counts[:, zero].sum(axis=1)
    ctrl0 = controls[zero].sum()
    if ctrl0 <= 0:
        return np.full(counts.shape, np.nan)
    expected = target0[:, None] / ctrl0 * controls[None, :]
    out = 1 - counts / (expected + 0.1)
    out[target0 <= 0] = np.nan
    out[:, controls <= 0] = np.nan
    return out


def audit(counts, cm, ko, times):
    zero = times == 0
    late = np.argsort(times)[-2:]
    pooled = normalized(counts, cm + ko, zero)
    cm_only = normalized(counts, cm, zero)
    ko_only = normalized(counts, ko, zero)
    late_pooled = pooled[:, late].mean(axis=1)
    shift = np.abs(cm_only[:, late].mean(axis=1) - ko_only[:, late].mean(axis=1))
    alternatives = []
    for index in np.where(zero)[0]:
        keep = zero.copy()
        keep[index] = False
        alternatives.append(normalized(counts, cm + ko, keep)[:, late].mean(axis=1))
    alt = np.array(alternatives)
    baseline_shift = np.max(alt, axis=0) - np.min(alt, axis=0)
    return {
        "f": pooled, "f_cm": cm_only, "f_ko": ko_only,
        "late": late_pooled, "control_shift": shift,
        "baseline_shift": baseline_shift,
        "late_pair_gap": np.abs(pooled[:, late[0]] - pooled[:, late[1]]),
        "baseline_reads": counts[:, zero].sum(axis=1),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path(__file__).parent)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "package")
    args = parser.parse_args()
    src, out = args.source, args.output
    out.mkdir(exist_ok=True, parents=True)
    z = zipfile.ZipFile(src / "count_data.zip")
    manifest = json.loads((src / "figshare.json").read_text())
    expected = {x["name"]: x["computed_md5"] for x in manifest["files"]}
    checksums = {}
    for name in ["count_data.zip", "HiTBiP_libraries.xlsx"]:
        digest = hashlib.md5((src / name).read_bytes()).hexdigest()
        assert digest == expected[name], (name, digest, expected[name])
        assert digest == SOURCE_MD5[name], (name, "Source has changed from the pinned analysis.")
        checksums[name] = digest
    for name, expected_hash in CASLAND_SHA256.items():
        digest = hashlib.sha256((src / "Casland" / name).read_bytes()).hexdigest()
        assert digest == expected_hash, (name, "Repository file differs from reviewed, pinned version.")
        checksums["Casland/" + name + " (sha256)"] = digest
    annotation = pd.read_excel(src / "HiTBiP_libraries.xlsx", sheet_name="subpool_array")
    # Original library reuses some element numbers for complex variants. Our
    # prespecified WT/single-substitution subset must be uniquely addressable.
    annotation = annotation[annotation.mutation.map(lambda s: s == "WT" or bool(SINGLE.fullmatch(s)))].copy()
    annotation = annotation.set_index(["subpool_name", "element_number"], verify_integrity=True).sort_index()
    metadata = pd.read_csv(z.open("count_data/association/5nM/metadata.tsv"), sep="\t")
    metadata["column"] = metadata["column"].str.strip()
    ddg = pd.read_csv(src / "Casland/S17_F7C_ddG_perturbation.tsv", sep="\t")
    ddg = ddg[ddg.assay == "Productive binding"].set_index("name")["ddG_perturbation"].to_dict()
    guides, summaries, rawrows = [], [], []
    files = [n for n in z.namelist() if n.startswith("count_data/association/5nM/count_table.")]
    files.sort(key=lambda n: int(re.search(r"S(\d+)\.tsv", n)[1]))
    for filename in files:
        sid = re.search(r"(S\d+)\.tsv", filename)[1]
        full = pd.read_csv(z.open(filename), sep="\t")
        full.columns = full.columns.str.strip()
        own = full[full.subpool_name == sid].copy()
        meta = metadata[metadata.subpool_name_sgrna == sid].sort_values(["time", "timepoint"])
        assert len(meta) in (15, 16) and sum(meta.time == 0) in (2, 3)
        columns = meta.column.to_list()
        counts = own[columns].to_numpy(dtype=float)
        assert np.isfinite(counts).all() and (counts >= 0).all()
        assert (counts == counts.astype(int)).all()
        cm = counts[own.mutation.to_numpy() == "CM"].sum(axis=0)
        ko = counts[own.mutation.to_numpy() == "KO"].sum(axis=0)
        assert len(own[own.mutation == "CM"]) == len(own[own.mutation == "KO"]) == 1
        times = meta.time.to_numpy(dtype=float)
        calc = audit(counts, cm, ko, times)
        # Publication includes separate negative controls. They audit normalization,
        # but are not candidate binders in the target review queue.
        candidates = ~own.mutation.isin(["CM", "KO"]).to_numpy()
        review = ((calc["baseline_reads"] < 100) | (calc["control_shift"] > .15)
                  | (calc["baseline_shift"] > .15) | (calc["late_pair_gap"] > .15)
                  | ~np.isfinite(calc["late"])) & candidates
        summary = {
            "id": sid, "name": str(annotation.loc[(sid, "E1"), "target_name"]),
            "all_targets": int(candidates.sum()), "timepoints": len(meta), "baseline_samples": int(sum(times == 0)), "low_baseline_100": int(((calc["baseline_reads"] < 100) & candidates).sum()),
            "control_sensitive_015": int(((calc["control_shift"] > .15) & candidates).sum()),
            "baseline_sensitive_015": int(((calc["baseline_shift"] > .15) & candidates).sum()),
            "late_pair_gap_015": int(((calc["late_pair_gap"] > .15) & candidates).sum()),
            "review_100_015": int(review.sum()),
            "cm_total": int(cm.sum()), "ko_total": int(ko.sum()),
        }
        summaries.append(summary)
        wt_seq = str(annotation.loc[(sid, "E1"), "cas_target"])
        guide = {"id": sid, "name": summary["name"], "sequence": wt_seq,
                 "times": times.tolist(), "timepoint_ids": meta.timepoint.astype(int).tolist(),
                 "controls": {"CM": cm.astype(int).tolist(), "KO": ko.astype(int).tolist()}, "targets": []}
        for i, (_, row) in enumerate(own.iterrows()):
            mut = row.mutation
            if mut != "WT" and not SINGLE.fullmatch(mut):
                continue
            seq = str(annotation.loc[(sid, row.element_number), "cas_target"])
            model_name, model_value = "perfect_target", 0.0
            if mut != "WT":
                pos, alt = mut.split(":")
                base = wt_seq[int(pos) + 20]
                classes = [("sub_complement", COMPLEMENT), ("sub_transition", TRANSITION), ("sub_opposing", OPPOSING)]
                label = next(label for label, trans in classes if base.translate(trans) == alt)
                model_name = f"{label}|1|{pos}"
                model_value = ddg.get(model_name)
                assert model_value is not None, model_name
            target = {"id": row.element_number, "mutation": mut, "sequence": seq,
                      "counts": counts[i].astype(int).tolist(), "f_bound": floats(calc["f"][i]),
                      "f_cm": floats(calc["f_cm"][i]), "f_ko": floats(calc["f_ko"][i]),
                      "late": number(calc["late"][i]), "baseline_reads": int(calc["baseline_reads"][i]),
                      "control_shift": number(calc["control_shift"][i]),
                      "baseline_shift": number(calc["baseline_shift"][i]),
                      "late_pair_gap": number(calc["late_pair_gap"][i]),
                      "model_name": model_name, "model_ddg": number(model_value)}
            guide["targets"].append(target)
            for j in range(len(meta)):
                rawrows.append([sid, guide["name"], row.element_number, mut, seq, int(meta.timepoint.iloc[j]), times[j], int(counts[i, j]), int(cm[j]), int(ko[j])])
        assert len(guide["targets"]) == 61, (sid, len(guide["targets"]))
        guides.append(guide)
        print(sid, len(own), summary["review_100_015"], flush=True)
    assert len(guides) == 91
    flat = [t for g in guides for t in g["targets"]]
    def flags(t, floor=100, tolerance=.15):
        return t["baseline_reads"] < floor or any(t[x] is None or t[x] > tolerance for x in ["control_shift", "baseline_shift", "late_pair_gap"])
    coverage = {"libraries": len(guides), "drilldown_targets": len(flat), "raw_count_measurements": len(rawrows),
                "all_non_normalizer_targets": sum(s["all_targets"] for s in summaries),
                "all_review_100_015": sum(s["review_100_015"] for s in summaries),
                "drilldown_review_100_015": sum(flags(t) for t in flat),
                "drilldown_review_30_030": sum(flags(t, 30, .30) for t in flat),
                "drilldown_review_300_010": sum(flags(t, 300, .10) for t in flat)}
    data = {"meta": {"title": "Filter-binding control and baseline sensitivity audit", "assay": "dCas9 association", "concentration_nM": 5,
                      "coverage": coverage, "default_baseline_floor": 100, "default_tolerance": .15,
                      "units": {"f_bound": "fraction depleted (not fitted occupancy)", "model_ddg": "kT, productive-binding energy penalty"},
                      "subset_rule": "All 91 sublibraries; drilldown includes WT and all 60 single substitutions at protospacer positions -20 through -1. All other targets are summarized, not in the drilldown.",
                      "normalization": "1 - target_count / (sum(target_zero)/sum(control_zero)*control_count + 0.1); metadata-listed zero-time points pooled (3, except S32 has 2); control=CM+KO",
                      "threshold_note": "Exploratory triage settings, not the authors' published QC and not a pass/fail assay acceptance criterion.",
                      "repository": "https://github.com/augustboyle/Casland", "commit": COMMIT,
                      "data_doi": "10.6084/m9.figshare.12526157.v2", "data_license": "CC BY 4.0", "code_license": "MIT"},
            "summaries": summaries, "guides": guides}
    (out / "data.json").write_text(json.dumps(data, separators=(",", ":"), allow_nan=False))
    with (out / "drilldown_counts.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sublibrary", "guide_name", "element", "mutation", "target_sequence", "timepoint_id", "time_min", "count", "CM_count", "KO_count"])
        writer.writerows(rawrows)
    with (out / "library_audit.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader(); writer.writerows(summaries)
    (out / "source_checksums.json").write_text(json.dumps(checksums, indent=2))
    (out / "findings.json").write_text(json.dumps(coverage, indent=2))
    print(json.dumps(coverage, indent=2))


if __name__ == "__main__":
    main()
