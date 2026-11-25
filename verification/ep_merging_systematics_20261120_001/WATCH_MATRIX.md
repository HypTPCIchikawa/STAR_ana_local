# Aggregate production completion and comparison watcher

`watch_production_matrix.py` watches the seven new case submissions as one study. It does not submit, cancel, resubmit, merge, remix or retry an analysis. It uses the `parse_events` function in `script/watch_job_notify.py`; that file is loaded without creating a Python cache file and its SHA256 is pinned. Queue state is obtained in one `condor_q` query for the seven explicit SUMS IDs per poll.

The parent launcher supplies `production_matrix.json` after all seven launch states have reached `submitted_verified` with `actualJobs:20`:

```json
{
  "schemaVersion": 1,
  "cases": [
    {
      "caseId": "ep1_mergingOff",
      "launchDirectory": "production/ep1_mergingOff",
      "config": "production/ep1_mergingOff/config.json"
    }
  ]
}
```

The real matrix must contain exactly the seven IDs in `cases.json:productionCases`, with no baseline row. Relative paths resolve against the matrix directory. `launch_state.json` must contain `phase`, `jobid`, `actualJobs`, `runmeta`, `eventLog`; the existing launch-config schema supplies `inputList`, `inputListSha256`, `sources`, `anaName`, `mainconf`, `outputDirectory`, `outputFileStem` and `submitHost`. Submitted XML snapshot SHA/size, native 20 input lists and their exact source order are verified independently. Existing EP8/cut-Off baseline state/config and completed CF input-stat proof are explicit references in `watch_matrix_settings.json`.

A read-only check of the prior actual baseline run is available immediately:

Before any submission, `--validate-settings --matrix PATH --settings PATH` validates settings, matrix/config references, configured file fingerprints and the 20 source stats. It does not read launch states, query Condor, inspect outputs, run postprocessing or send mail. The full matrix declares `projectRoot`, `studyDirectory`, `casesManifest` and `casesManifestSha256`. The `additionalPinFiles` setting explicitly includes the parent preparation and matrix-runner helpers. Freeze matrix, launch states, configs and helper files once normal watching starts.

```bash
python3 rootfile/femto_lambda_ep_merging_systematics_20261120_001/watch_production_matrix.py \
  --settings rootfile/femto_lambda_ep_merging_systematics_20261120_001/watch_matrix_settings.json \
  --check --check-baseline
```

A matrix check uses `--check --matrix PATH`. Check mode performs one event/queue pass. If all jobs are terminal, it also verifies logs, native input lists and ROOT headers/stats. It prints JSON and creates no files, mails, plots or analysis processes. Check-mode errors return nonzero and send no notification.

Start the durable watcher from the project directory on the submit host, using an explicit Python path and a new state directory. The following is a command template; the parent launcher starts it after the real matrix is complete:

```bash
nohup /cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/python-3.8.13-myixjs7jesskpbhwiowibizezfh52rky/bin/python3 -u \
  rootfile/femto_lambda_ep_merging_systematics_20261120_001/watch_production_matrix.py \
  --settings rootfile/femto_lambda_ep_merging_systematics_20261120_001/watch_matrix_settings.json \
  --matrix rootfile/femto_lambda_ep_merging_systematics_20261120_001/production_matrix.json \
  --state-dir rootfile/femto_lambda_ep_merging_systematics_20261120_001/aggregate_watch_001 \
  > rootfile/femto_lambda_ep_merging_systematics_20261120_001/aggregate_watch_001.log 2>&1 < /dev/null &
```

The watcher polls every 300 seconds. Each new case needs exactly 20 submitted event identities, all terminal with `exit_0`, and an empty matching queue. A held worker, removal, nonzero/unknown exit, input drift or a verification/command failure records `attention.json` and stops. After all 140 new workers become terminal, the watcher waits 300 seconds and confirms the condition again before output verification. This settle interval accommodates copy/log visibility; it does not retry an analysis.

For each of the 140 new and 20 reused baseline outputs, it requires:

- Exact source period/order from the native input list and exactly one ROOT for each index 0–19 in the exact SUMS filename stem.
- Exactly one replay marker with job ID and `requested=completed=available=expectedPeriodEvents`, `status=0`.
- Maker `status=completed`, exact reconstructed/input event counts, chain nStOK equal to expected events and zero nStWarn/EOF/Err/Fatal.
- Normal native exit/END, no configured fatal/ROOT-I/O signature in stdout or stderr; empty stderr is valid.
- A nonempty ROOT with `root` header and stable byte size/inode/device/mtime. ROOT source-ID/settings metadata is then verified by the histogram-only CF postprocessor.

All cases require 1,819,447,887 completed stored events. Original first-stage period stats, reused baseline output stats, watcher/postprocess source SHA256, settings/cases/matrix, launch states/configs/runmeta/submitted XML and native lists are pinned. Inputs are checked again before and after postprocessing. Logs are read once at final verification, with SHA256 computed during that read. There is no TTree entry loop.

Under the state directory, `completion_verified.json` contains the 160-worker proof; `output_lists/*.list` and `case_outputs.json` provide the strict postprocessor input schema. `postprocess_attempt.json` is persisted exclusively before launching the single comparison command. The comparison outputs go to `cf_comparison/`, followed by `postprocess_verified.json`. The user-facing figure contract and command are in `CF_POSTPROCESS.md`.

Successful notification additionally requires exactly one inclusive 0–60% CF with positive SE/ME normalization counts in each of the eight cases. Both postprocessor and watcher enforce the `requiredValidCentralities` setting. Empty peripheral or individual-plane curves remain explicitly unavailable.

Only one minimal Japanese completion or attention email is attempted, to the explicit address in the settings. There is no per-case mail. The message contains no job IDs, paths, counts or technical diagnostics. The watcher uses a process lock and exclusive `delivery_attempt.json` persisted before `/usr/sbin/sendmail -i -t`; `sent.json` records local MTA acceptance. A failed/uncertain send is not automatically resent. Existing attention, delivery-attempt or sent states block automatic restart. Failures remain in local JSON/logs, and the failure message does not claim that all jobs or physics validation completed.

Restart while only polling may reuse an unchanged startup pin record. Once postprocessing or mail has been attempted, any partial result requires manual inspection; the tool does not delete attempt gates or rerun incomplete postprocessing. `sent.json` proves local sendmail acceptance, not receipt in the destination mailbox.
