# Single-He4 laboratory-rapidity postprocessing and completion watch

These helpers operate only on saved pair histograms and small metadata. They do not read event-tree entries, rerun mixing, submit jobs or hadd outputs. The original EP study helpers are imported/included read-only and SHA-pinned for reproducibility.

`he4_y_postprocess_settings.json` fixes EP8, cut Off, depth100, continuous within each period, lab Y8 in `[-1,0)`, coarse Y4 by adjacent-two raw-count/Sumw2 sums, and H=`[0.4,0.5)` GeV/c. Single-particle He4 rapidity comes from its physical P4; it is distinct from pair rapidity. The pipeline retains and reports Y underflow/overflow. Four independent controls are normalized separately when possible: all Y including flow, sum of the eight in-range bins, Y<-1 and Y>=0. An empty H means unavailable/NA, never a fabricated zero CF or a fallback normalization.

Production output list must contain exactly twenty distinct ROOT paths in source-period order. Each output is matched to the corresponding original `periodMix100` result. Full/signal/leftSB/rightSB canonical SE/ME counts/errors/Entries are checked against that reference, while all new TH3 Y bins **including Y flow** must project back to their matching inclusive histogram. k* and centrality axes, explicit Y metadata, common event/EP-valid population and replay/mixing settings are checked before plotting.

```bash
python3 rootfile/lambda_he4_ylab_split_20261121_001/postprocess_he4_y.py \
  --settings rootfile/lambda_he4_ylab_split_20261121_001/he4_y_postprocess_settings.json \
  --outputs OUTPUTS.list --output-dir NEW_RESULT_DIRECTORY
```

For the same-prefix 400k pilot, add all four explicit arguments:

```text
--pilot --pilot-events 400000 --pilot-source FIRST_PERIOD.root --pilot-reference rootfile/femto_lambda_period_mix100_20261119_001/pilot_B_period100.root
```

Every output directory must be new. Main files are `he4_y4_CF.pdf`, `he4_y8_CF.pdf`, `he4_y4_SE_ME_scaled.pdf`, `he4_y8_SE_ME_scaled.pdf`, `he4_y_flow_controls.pdf`, companion `cf_histograms.root`, `he4_y_cf_metrics.json`, input provenance and `completed.json`. CF errors hold alpha fixed, omit normalization/shared-event/shared-ME covariance, and are not a significance test. Native k* flow is listed separately and excluded from H and finite high-k windows. Windows outside the saved axis are unavailable.

## Watcher interface

Before submission, after `launch_config.json` exists:

```bash
python3 rootfile/lambda_he4_ylab_split_20261121_001/watch_he4_y.py \
  --launch-config rootfile/lambda_he4_ylab_split_20261121_001/launch_config.json \
  --launch-state rootfile/lambda_he4_ylab_split_20261121_001/launch_state.json \
  --settings rootfile/lambda_he4_ylab_split_20261121_001/he4_y_watch_settings.json \
  --validate-settings
```

This mode reads neither launch state, scheduler nor production outputs and sends no mail. After the one-shot submission is verified, use the same arguments with `--state-dir rootfile/lambda_he4_ylab_split_20261121_001/completion_watch_001` instead of `--validate-settings`. `--check` is a one-shot, read-only scheduler/completion check without plotting or mail.

The state may be `submitted_verified`, `monitor_starting` or `monitor_started`; owner phase updates do not alter the captured job identity. Config, runmeta, submitted XML, helper files, postprocess settings and twenty source ROOT stats are pinned. All twenty terminal exit0, queue empty, settling, exact native input order, log completion/counters, normal exits and ROOT headers must pass before postprocessing starts once. Comparison completion and unchanged output stats are required before the single minimal completion email. Errors create attention state and one minimal failure email. Reprocessing, resubmission, hadd and automatic retries/resends are disabled.

The parent launcher must keep the launch config and helper/configuration files unchanged while this watcher runs. A previous mail attempt or attention record requires inspection rather than an automatic restart. No watcher has been started by the preparation task.
