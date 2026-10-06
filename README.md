# Filter-binding control and baseline sensitivity audit

[Open the interactive demo](https://liuzhitong330.github.io/greenleaf-lab-demo/). Prepared by Cathy Liu with AI-assisted implementation and independent computational checks. This is a personal exploratory contribution, not an endorsed lab tool or evidence of hands-on filter-binding experience.

This independent analysis starts from deposited sequencing-count tables, not FASTQ files. It is a sensitivity audit, **not** a reproduction of the publication's likelihood fits, classification, fitted final occupancy, kinetic rates, or confidence intervals.

## Sources

- Boyle, Becker, Bai et al. (2021), *Science Advances* 7:eabe5496, DOI [10.1126/sciadv.abe5496](https://doi.org/10.1126/sciadv.abe5496).
- Author-deposited data, Figshare [10.6084/m9.figshare.12526157.v2](https://doi.org/10.6084/m9.figshare.12526157.v2), CC BY 4.0. Files: `count_data.zip` and `HiTBiP_libraries.xlsx`. All source authors retain attribution; this is a transformed analytical subset.
- Evan Boyle's [Casland](https://github.com/augustboyle/Casland), commit `9193bf679c2bab346d8c0f30a1f92e52c29c9bed`, MIT. The README explicitly associates the repository with this paper. Its productive-binding energy-perturbation table is used for 60 protospacer single-mismatch classes. The reviewed original command-line script was actually run for the L1 guide; all 60 generated sequences and energies were checked against the derived data. Other guide mappings apply the same substitution rules.

## Coverage

All 91 deposited 5 nM association sublibraries are audited. The detailed explorer is a **prespecified subset**: one perfect target plus all 60 single substitutions at positions -20 through -1 per sublibrary, yielding 5,551 target records and 88,755 count observations. It excludes context variants, PAM variants, bulges, and multi-mismatches from drilldown. Summary CSV includes all 54,167 own-sublibrary target rows except the CM and KO normalization controls. These are sequence-library measurements, not 54,167 biological replicates. Different sublibraries and time points are not independent replicates.

Only count-table rows whose `subpool_name` equals the assayed sublibrary are selected. The tables also contain rows for other sublibraries; those are not incorrectly treated as additional binding targets. Original annotations reuse some element numbers for complex variants, so the WT/single-substitution subset is verified to have a unique `(subpool,element)` identity.

The deposited metadata, not nominal methods timings, determines sample times and order. S32 has 15 listed samples and two zero-time controls; the other 90 sublibraries have 16 samples and three zero-time controls. An unlisted S32 sample is not imputed. The reason for its metadata omission is unknown.

## Calculation

For each target and normalization-control choice C (CM+KO pooled, CM alone, or KO alone), let `target0=sum(target counts at listed zero-time samples)` and `C0=sum(C counts at those samples)`. At time t:

```
expected_count(t) = target0 / C0 * C(t)
fraction_depleted(t) = 1 - observed_count(t) / (expected_count(t) + 0.1)
```

The `+0.1` count pseudocount follows the paper's plotted fraction-bound equation. Pooling available zero-time counts is an explicit choice in this audit, not a claim about the authors' unpublished baseline implementation. A zero target baseline, zero control baseline, or zero time-point control count is unidentifiable and retained as null. Negative estimates and estimates above one are not clipped. This is a count-normalized depletion proxy for binding; it is not a fitted equilibrium occupancy. No kinetic model is fit here.

Late depletion is the mean of the final two chronological time points (59 and 60 minutes). Three transparent sensitivity metrics are reported in **fraction units**:

1. **Control sensitivity:** absolute difference between the late CM-only and KO-only estimates. It probes dependence on normalization choice; it cannot identify which control is correct.
2. **Baseline sensitivity:** range of late pooled-control estimates when each available zero-time sample is omitted in turn. It measures baseline dependence, not a confidence interval or biological variability.
3. **Late-pair discrepancy:** absolute difference between the two latest pooled-control estimates. It may reflect counting noise, baseline/control noise, timing or biology. It is not evidence of a specific equipment failure.

A review flag means baseline reads below the user-selected floor, any sensitivity metric above the selected tolerance, or a missing sensitivity estimate. The default is 100 pooled baseline reads and 0.15 (15 percentage points). These are **exploratory settings**, not published QC cutoffs, validated acceptance criteria, proof that a target failed, or proof that unflagged targets are correct. Strictly `>` tolerance and strictly `<` baseline floor are used. Rounded display data are sufficient to reproduce the selected review counts. The full-source generator computes before rounding.

## Repository-derived model annotation

Casland supplies relative perturbations to productive-binding energy barriers in kT. Values are joined by protospacer position and substitution class (complement, transition, opposing base), with the original 23-base non-template target orientation and PAM retained. The perfect target is the zero-perturbation reference. These numbers are **not** binding probabilities, 5 nM occupancies, independent validation, or a predictor of instrument performance. The source model was fit using related source experiments. They can help select a prospective challenge panel spanning mismatch types, followed by independent validation.

## Reproduce

Use Python 3.11+ with `numpy`, `pandas`, and `openpyxl`; tested package versions are in `analysis/requirements.txt`. From this repository root, download original sources into an empty directory and clone the pinned author repository. No private data or credentials are needed.

```sh
python -m pip install -r analysis/requirements.txt
mkdir source
curl -L --fail https://api.figshare.com/v2/articles/12526157/versions/2 -o source/figshare.json
curl -L --fail https://ndownloader.figshare.com/files/23321960 -o source/count_data.zip
curl -L --fail https://ndownloader.figshare.com/files/24177854 -o source/HiTBiP_libraries.xlsx
git clone https://github.com/augustboyle/Casland.git source/Casland
git -C source/Casland checkout 9193bf679c2bab346d8c0f30a1f92e52c29c9bed
python analysis/build_data.py --source source --output reproduced
python analysis/test_data.py --package reproduced --casland source/Casland
```

`build_data.py` verifies the source archive/workbook MD5s against the Figshare manifest. `source_checksums.json` records the independently downloaded digests. `test_data.py` independently recalculates late depletion, control sensitivity, baseline omission ranges, and review flags using scalar arithmetic, checks all 88,755 raw subset CSV rows, and invokes the reviewed author script for an exact sequence/energy cross-check. It does not execute unreviewed downloaded code.

Outputs: `data.json` (browser data), `drilldown_counts.csv` (auditable raw-count subset), `library_audit.csv` (all-target library summary), `findings.json` (coverage and settings-dependent counts), `validation.json` (checks), and `casland_L1_binding_cli.tsv` (actual upstream CLI output).

Published outputs are in `data/`; the count subset is distributed as `drilldown_counts.csv.gz` to keep downloads small. The test script accepts this compressed form. Serve the static page with `python -m http.server 8000` from the repository root; opening `index.html` directly with `file://` may block its data fetch.

## Downloads and prospective use

- [Auditable raw-count subset](data/drilldown_counts.csv.gz)
- [Whole-library summary](data/library_audit.csv)
- [Computational validation record](data/validation.json)
- [Blank paired-run log](run-log-template.csv)

The website can export its current target queue, thresholds and review reasons as CSV. The blank run log contains no fabricated experimental records. Equipment validation remains a proposed experiment requiring independent days, suitable controls, and prespecified acceptance margins.

Source data remain CC BY4.0 (Boyle etal.,2021); source Casland code/parameters retain the MIT attribution in `CASLAND_LICENSE`. This repository does not license or redistribute private applicant material.

## Reuse and limits

A next-run input needs target identity and sequence, sublibrary, protein concentration, sample/time-point IDs, actual times, target counts, and separate CM/KO control counts. New instrument experiments additionally need run/batch identity, well position, handling order, pressure/vacuum trace, coating age, recovery controls, and independent repeats. Those logs are absent from this analysis. A reusable next step is to rerun the same normalization sensitivities on a independently collected reference panel and compare review flags with experimental outcomes. The current data do not establish an apparatus defect, a repair, or a validated performance improvement.
