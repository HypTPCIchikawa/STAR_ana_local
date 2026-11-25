# Period-continuous ME with depth100

[Task/status](../../Todo/2026/11/19/period_mix100_production.md).

The submitted campaign is SUMS `2D07F680577123CABC16FDD9159CFD9E`: one
validated period ROOT per worker, 20 workers. Both pool depth limits are 100,
producer segment boundaries retain the pool, and cent9/vertex/EP binning and
bidirectional buffer mixing are preserved. The 34 unrecovered producers and
81 uncalibrated runs remain deferred.

- `evidence/pilot_comparison.json`: 400k-event real-tree comparison of continuous100,
  reset100 and continuous20. All 607 SE histograms match exactly; He4 signal ME
  entries are 116183, 109474 and 24000. Three producer segments are traversed.
- `evidence/pilot_B_ep8_plane_sum.json`: 16 inclusive/EP-sum comparisons, all exact.
- `evidence/dictionary_fresh_process_smoke.log` and the before/after LinkDef snapshot:
  fresh normal replay succeeds after explicitly registering the primitive-vector
  dictionaries. Failed initial trial evidence remains in the parent study's
  `initial_dictionary_failure/`; it is not a successful pilot.
- `snapshot/`: pinned textual source/config/helper snapshots. Binary libraries
  and calibration ROOTs are recorded by SHA in the launch config and actual
  submit-time ZIP/runtime archives; they are not duplicated here.
- `submitted/`: actual runmeta/configlog/joblistlog from this submission.
- `preserved_manifest.json`: copied-file SHA256 and gzip flags, verified by readback.

The parent study is `rootfile/femto_lambda_period_mix100_20261119_001/` and the
outputs are `rootfile/auau3p85_anaFemtoLambda_4He_tree_periodMix100/production20261119/`.
All 20 production jobs and output metadata checks pass. The
`completion_20261120/` directory preserves the full completion-log audit,
bounded ROOT check, frozen output list and completion-mail receipt. All
1,819,447,887 stored events are accounted for; output ROOTs total 11,491,533,504
bytes. The ROOT audit reads only 40 scalar-ledger endpoint rows and does not
perform a full event replay. The original `evidence/` launch snapshots remain
historical. No automatic resubmit or hadd is enabled.
