# Raw-signal CF histogram postprocessing

`postprocess_cf.py` reads only four saved signal pair histograms per ROOT and small provenance objects. It does not read audit/rich-tree entries, remix events, merge ROOT files or submit jobs.

The centrality contract follows `star_local/analysisnote/20261112/ep_and_cf_combination_review.md` and the existing `checkHistAnaFemtoLambda.C` / EP-rapidity CF study: add raw SE and ME counts within the selected cent9 range, then normalize once with `alpha=SE(H)/ME(H)`, `CF=SE/(alpha*ME)`, `H=[0.4,0.5)` GeV/c. CFs are not averaged. The signal remains before sideband subtraction/purity correction. Plane-inclusive CF uses the sum of raw plane counts, not an average of plane CFs. Each diagnostic plane gets its own normalization. All accepted cent9 0–8, 0–60%, 0–20%, 0–10%, 10–20% and each independent cent9 bin are in the settings file.

Run from the project directory. Both commands require a **new** output directory.

```bash
python3 rootfile/femto_lambda_ep_merging_systematics_20261120_001/postprocess_cf.py \
  --cases rootfile/femto_lambda_ep_merging_systematics_20261120_001/cases.json \
  --settings rootfile/femto_lambda_ep_merging_systematics_20261120_001/cf_postprocess_settings.json \
  --mode baseline \
  --output-dir rootfile/femto_lambda_ep_merging_systematics_20261120_001/cf_baseline_ep8Off_NEW
```

Baseline mode uses only the explicitly pinned existing EP8/cut-Off 20 output paths in `cases.json`. It does not create panels for absent cases. The comparison mode requires an output-list mapping for **all eight** cases:

```json
{
  "schemaVersion": 1,
  "cases": {
    "ep1_mergingOff": {"outputList": "output_lists/ep1_mergingOff.list"},
    "ep1_mergingOn": {"outputList": "output_lists/ep1_mergingOn.list"},
    "ep8_mergingOff": {"outputList": "output_lists/ep8_mergingOff.list"},
    "ep8_mergingOn": {"outputList": "output_lists/ep8_mergingOn.list"},
    "ep16_mergingOff": {"outputList": "output_lists/ep16_mergingOff.list"},
    "ep16_mergingOn": {"outputList": "output_lists/ep16_mergingOn.list"},
    "ep32_mergingOff": {"outputList": "output_lists/ep32_mergingOff.list"},
    "ep32_mergingOn": {"outputList": "output_lists/ep32_mergingOn.list"}
  }
}
```

Relative list paths resolve against the mapping JSON directory. Relative ROOT paths resolve against their list directory. Each list must contain exactly 20 distinct regular ROOT paths, one path per line, ordered like `periods.list` / `input_manifest.json`. The baseline list must exactly match `existingBaselineOutputs`. Cross-case duplicate paths/inodes are rejected. An arbitrary ROOT directory is never globbed to construct coverage.

```bash
python3 rootfile/femto_lambda_ep_merging_systematics_20261120_001/postprocess_cf.py \
  --cases rootfile/femto_lambda_ep_merging_systematics_20261120_001/cases.json \
  --settings rootfile/femto_lambda_ep_merging_systematics_20261120_001/cf_postprocess_settings.json \
  --outputs rootfile/femto_lambda_ep_merging_systematics_20261120_001/case_outputs.json \
  --mode compare \
  --output-dir rootfile/femto_lambda_ep_merging_systematics_20261120_001/cf_eight_cases_001
```

Validation requires completed replay status, one correct source period per output, exact expected stored-event processing, equal reconstructed selection/stored-EP calibration/reduced settings, common pool100 settings except the declared EP count, false/false segment-reset/original-mixing, and equal EP-valid event counters per period. New cases require explicit cut-enabled and pair-EP-count metadata. Only the pinned original EP8/cut-Off SUMS result may omit those two new metadata objects. Summing the regular EP bins must reproduce the inclusive signal histogram, including k* flow, within the explicit count tolerance. Nonzero centrality/plane flow fails. Input source ROOT stats and config/helper/list SHA256 are checked again after plotting.

Outputs:

Each selected case requires a valid inclusive 0–60% normalization with positive SE/ME counts (`requiredValidCentralities`). Empty peripheral/plane curves remain unavailable. New pair-cut definitions must agree on proton–He4 target, strict AND thresholds 0.05/0.05, r=1.4 m, same SE/ME acceptance, unreachable-retained policy and unchanged He3/pion policy. Inclusive SE counts must agree across EP choices at fixed cut state. The verified historical baseline is `cf_baseline_ep8Off_002/`; it pins its calculation's helper version. Later preflight/main-curve guards do not change its recorded counts or CF values.

| File | Content |
|---|---|
| `cf_comparison.pdf` | Baseline-only cent panels, or 4 EP panels per centrality with cut Off/On overlays |
| `cf_plane_diagnostic.pdf` | Each case's separately normalized plane CF, with inclusive reference; default 0–60% |
| `raw_se_me_scaled.pdf` | SE/SE(H) and ME/ME(H) shape overlays; raw retention is in JSON |
| `cf_comparison.json` | Per-curve alpha, SE/ME norm and full-axis counts, k* flow, display bins, high-k window ratios/means; On/Off raw retention |
| `cf_histograms.root` | Postprocessed raw signal SE/ME and CF graphs; invalid CF points are absent, not zero-filled |
| `input_provenance.json`, `raw_histograms.json`, TSV/logs | Exact input stats/config pins, raw histogram bins, validation and plot inputs |
| `completed.json` | Success marker, unchanged-input checks and output hashes; `failed.json` marks failure |

High-k `integralCF` uses the ratio of window-integrated raw counts with the curve's alpha. `arithmeticMeanOfValidDisplayBins` is a separate descriptive mean. A 1.0–1.5 GeV/c request outside a 0–1 native axis is explicitly unavailable; the tool does not clip or extrapolate. Native k* under/overflow is reported separately and is excluded from normalization/CF. Raw scaled plots are shape comparisons, not efficiency/retention measurements.

Display errors use histogram Sumw2 propagation while alpha is fixed. They omit normalization covariance, shared-event/shared-ME correlations and cross-case covariance. Window comparisons are descriptive; the helper performs no significance test. It does not infer new mixing results from existing histograms. EP1 is required to share the same stored EP-valid sample.
