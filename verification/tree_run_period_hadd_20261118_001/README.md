# Tree run-period hadd evidence

- [Task and live-state locations](../../Todo/2026/11/18/tree_run_period_hadd.md).
- [Grouping report](../../analysisnote/20261118/tree_run_period_grouping_20261118.md).
- [Missing-input report](../../analysisnote/20261118/tree_missing_inputs_20261118.md).
- `execution/`: frozen 4024-source snapshot/list, operation config, pinned merge
  driver/ROOT helpers, two-source merge/check proof and quota/launch records.
  `execution/preserved_manifest.json` records the SHA256 of the uncompressed
  original bytes. Artifacts larger than 1 MiB are saved with gzip.
- `grouping/`: 10/20/30 comparison, adopted 20-period table and snapshot proof.
- `failure_coverage.json`, `failed_worker_inputs.tsv`, `run_coverage.tsv`, and
  `current_prior_comparison.tsv`: final 34-producer loss audit.

Production ROOTs, merged ROOTs and bulk logs remain in the parent study paths.
These copies describe the launch and bounded IO validation; they do not certify
completion of the asynchronous 20-group merge. Check live `merge_status.json`
and final `completed.json` in the parent study for completion evidence.
