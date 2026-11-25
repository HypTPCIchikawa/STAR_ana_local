# Run-period grouped reduced trees

This study concatenates completed `anaFemtoLambda_4He_tree` producer files into
20 contiguous run periods. `merge_config.json` fixes the source snapshot, grouping,
ROOT helper hashes, runtime, resource settings and output directory. The production
sources, native SUMS lists and runtime/calibration records are retained.

## Inputs and grouping

- SUMS: `9037A6F05FAD08FC8E9BD9286C8A38A5`, 4,058 submitted producers.
- `source_snapshot_final.json`: 4,024 completed producer ROOTs, including 2869.
- Successful input: 72,442 PicoDst files, 387 represented runs, 1,842,806,823
  input entries and 1,819,447,887 saved events, from producer completion logs.
- `grouping_analysis/computed_final_run_cohorts/groups_20/assignment.json` and
  `summary.csv`: whole producers in minimum-run cohorts, 185–240 inputs per
  period, 22.805–24.236 GB source bytes per period.
- The grouping uses PicoDst run IDs rather than generated ROOT modification times.
  Cross-run producers produce small overlaps in 18 runs at adjacent boundaries;
  no event tree or producer segment is split. The minority side contains 175 Pico
  files, 0.2416% of the successful 72,442 Pico files.
- `source_snapshot.json` retains the earlier 4,023-file snapshot before producer
  2869 completed; it is historical, not the operational merge input.

## Running and checking

The tmux session is `he4_tree_period_hadd_20261118_001` on `starsub03`. The driver
is `run_grouped_merge.py --config merge_config.json` with absolute paths in the
launch record. `tmux_launch.json` stores the exact command. Do not start another
instance: the driver has a process lock and refuses existing status/output paths.

`merge_status.json` and `tmux_merge.log` show current progress. For each group,
`operations/group_XX/` contains the frozen inputs, source catalog, independent
snapshot match, quota check, merge report, output check and publication record.
The source stat fingerprint is checked before and after work. Producer segment
keys, branch/QA schemas and provenance are checked across groups.

ROOT 5.34 / SL24y runs through the production mainconf and batch wrapper. The
project `vector<TVector3>` dictionary is loaded. The fast merger keeps compression
6 and disables local input copies. The same ROOT process explicitly sets maximum
tree size to 1 TB. All group inputs fit the open-file allowance; one merge runs at
a time. The quota guard requires 110% of remaining source bytes plus 50 GB below
the soft limit after subtracting in-doubt usage.

Outputs go to `rootfile/auau3p85_anaFemtoLambda_4He_tree/run_periods20261118_001/`.
A `.partial.root` becomes a formal output only after independent checks. The
driver publishes without replacing an existing destination. Only after every
group succeeds does it write `grouped_roots.list`, `completed.json` and the global
source catalog validation. An error leaves `attention_required` and retained
partial evidence for inspection; no input is skipped automatically.

Checks cover all tree schemas/counts, ordered segment metadata, all 108 QA
histograms including under/overflow, errors and Entries, three identity samples
per segment and bounded official Reader payload samples. QA float allowance is
0 ULP, established by the two-source real merge in `validation/`. No complete
physics replay or exhaustive event-payload corruption scan is claimed.

## Missing inputs and second stage

`failure_audit/failure_coverage.json` and its TSV tables compare the 34 missing
producers with previous EP1/EP8 inputs. 33 identical offending Picos are identified;
one repeated Make failure lacks a specific offending Pico. The missing producers
contain 612 of 73,054 Picos (0.837736%). Known Entries total at least 8,594,396;
302 Picos have unknown Entries. The known missing-input fraction is at least
0.464210%, with substantial run concentration. This is not an accepted-event or
pair-loss fraction. Run 22162008 has only a measured zero-entry Pico.

Stage 2 is not launched by this driver. Standard replay currently resets the
pool at each original producer segment. Period-wide mixing requires a reviewed
variant with `resetPoolAtSegmentBoundary:false`, `requireOriginalMixing:false`,
and one period per stage-2 job. Hadd alone does not change the pool-reset policy.
