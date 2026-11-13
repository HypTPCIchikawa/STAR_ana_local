# Quota review for the 388-run EP comparison

Status: projection and same-input 10000-event pilot comparison complete. A packed 25-file production worker has not been measured. No job submitted and no file deleted.

The current user soft-quota headroom after inDoubt is **286.056 GB**. GB means 10^9 bytes. User usage is 2.302 TB, soft 2.588 TB, hard 3.138 TB. File-count soft/hard limits are 0/0 (not imposed). Group rhstar has no block/file limit; star-pwg fileset headroom is about 214 TB and is not the binding constraint.

Prior hist-only campaign: 26,945 successful worker ROOTs, 160.245 GB combined, 5.947 MB mean, 11.019 GB final hadd. The old individual shards are now absent; the manifest is the immutable size evidence.

## Two-campaign projection (EP1 plus EP8)

| Files/job | Jobs/campaign | Reference effective-fixed model | All new arrays charged raw | Double-reference stress + raw new arrays |
|---:|---:|---:|---:|---:|
| 5 | 14611 | 303.7 GB | 420.1 GB | 601.9 GB |
| 10 | 7306 | 203.7 GB | 261.9 GB | 362.7 GB |
| 25 | 2923 | 143.6 GB | 166.9 GB | 219.2 GB |

Each number includes all source worker ROOTs, four retained merge-round equivalents (39.8 GB/campaign), four possibly nonempty 256-KiB small records per worker, a 256-KiB ROOT-allocation rounding reserve per worker, and 2 GB/campaign for runtime/runmeta/audit records. Final hadd is published by hardlink and is not counted twice.

The new-array raw estimate is 3.982896 MB/worker: 1.706496 MB for 8 new TH3D content arrays plus 2.276400 MB for the prior 32 YLab, 6 pion-merging and 3 EP additions. Actual ROOT compression is expected to reduce this but is not assumed in the raw column.

The decomposition F=(R−M)/(N−1), E=M−F gives effective F=5.538370 MB/old worker and 9.954 GB for E in the covered inputs. It is an empirical model; histogram occupancy/compression is nonlinear. The double-reference stress case doubles the entire reference ROOT component before adding the raw new arrays.

**Capacity recommendation: 25 files/job**, 2923 jobs/campaign for a flat list. Two campaigns project to 166.9 GB with every new content array charged uncompressed; even double-reference stress projects 219.2 GB, leaving about 66.9 GB to current soft quota. Actual packed-worker output is not yet measured. 5 files/job does not have enough two-campaign margin; 10 is less robust. If grouping is reset at run boundaries, up to 387 extra workers are possible and must be counted from the actual generated joblist.

## Runtime and scope

The preserved Condor submit sample has no explicit request_memory, request_disk or absolute walltime. Its PeriodicRemove removes CPU-inefficient jobs after 54000 sec (>15h and CPU fraction<.10), plus restart/hold conditions; it does not establish an absolute allowed runtime. See quota_condor_constraints.json. Worker input Pico files and runtime bundles are staged to worker SCRATCH, not retained as full copies in user GPFS. More files/job increases staging disk and retry scope and changes pool resets at job boundaries; EP1 and EP8 must use identical packing.

Own Condor queue on starsub03 had 0 jobs; a scoped own-process read found no root4star/hadd/submit/watcher processes at the check. Other submission hosts and independent filesystem writers are not guaranteed absent. Check quota again immediately before submission and retain a buffer.

No candidate store and no jointSparse are assumed. The existing candidate-store production around 2 TB is retained and has not been deleted to make room. No cleanup or ledger removal is part of this recommendation.

Machine-readable evidence: quota_user_raw.log, quota_user_command.json, quota_projection_models.json, quota_projection_pre_pilot.json, quota_condor_constraints.json.


## Same-input new-histogram pilot measurement

The existing EP8 10000-event ROOT is 4,039,238 logical bytes; the new pair-rapidity ROOT is 4,057,009 bytes, an increase of **17,771 bytes**. Both occupy **4,194,304 allocated bytes** on GPFS, so the measured allocation delta is zero for this pilot. The 8 new TH3D compressed TKey objects total 13,895 bytes (parent's dynamic audit). This is a small-statistics occupancy measurement, not an exact prediction for a 25-file worker. The conservative projection still charges all 3.983 MB of added content arrays per worker uncompressed.

The machine-readable result is `quota_pilot_storage.json`; all source/config/job settings remain owned by the parent implementation. A full-statistics submission is not part of this audit.
