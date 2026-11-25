# Compact reduced-tree grouping helpers

Run `run_grouped_tree_validation.C` through the existing batch-matched
`script/singularity_root_batch.sh` with the explicit tree-production mainconf.
The loader follows the repository reduced-tree IO test and loads the project
`vector<TVector3>` dictionary before compiling the study helper.

The nine explicit arguments are:

`(action, absoluteInputList, freshReport, stagedRoot, maxTreeBytes, expectedQa, sampleEvents, maxOpenedFiles, floatUlps)`

- `probe`: load dictionaries and report ROOT5's original/effective maximum tree
  size, effective open-file limit, and fast/local-copy configuration. Does not merge.
- `catalog`: check every source header, all producer metadata, all tree entry
  counts and branch schemas, and exactly `expectedQa` QA histogram headers.
  Producer segment keys must be unique; metadata counts must cover the complete
  event tree. Check each segment's first/middle/last identity and validate the
  first `sampleEvents` events with the official Reader. A report finishes only
  after every source succeeds. No complete event replay is performed.
- `merge`: merge complete source files in list order with `TFileMerger(false,true)`,
  fast method, unchanged compression, and a fresh `CREATE` output. The maximum
  tree size must exceed group source bytes plus 10% and 1 GB. Effective open-file
  capacity must cover the complete group plus two, avoiding hidden histogram
  merge batches. The helper never publishes, archives, deletes, or notifies.
- `check`: repeat the bounded source checks, require exact ordered producer
  metadata and all tree entry/schema sums in the staged output, replay all QA
  histograms with native ROOT Add, and compare every bin including under/overflow.
  QA Entries are exact. Float content tolerance is explicitly `floatUlps` ULPs
  relative to native Add (`0` means exact); the observed maximum and changed-bin
  count are reported. Double content and Sumw2 comparisons use an accumulation
  roundoff bound proportional to the number of inputs. Stored variances must be
  finite and nonnegative. This is IO validation, not an all-event physics replay.

Use a maximum such as 1 TB explicitly in the *same ROOT process* that merges.
Ordinary `hadd` in a different process does not inherit this setting. Keep inputs
whole; arbitrary event splitting can invalidate producer segment metadata.
Choose any nonzero Float ULP allowance from a reviewed small merge probe and
record it in the operation configuration; it is not a blanket relative tolerance.

After a full source catalog succeeds, `check_source_catalog.py` matches every
path, event counter and current filesystem fingerprint against `source_snapshot.json`.
It also checks global segment-key uniqueness and common branch/QA schemas.
The parent driver must verify input fingerprints before and after operations,
validate each staged output, and only then publish without clobbering.

No helper claims to detect an unsampled corruption inside an event basket or a
duplicate physical event outside the bounded Reader samples. Fast basket copying,
exact entry and segment preservation, and the sampled payload decoding form the
stated validation scope.
