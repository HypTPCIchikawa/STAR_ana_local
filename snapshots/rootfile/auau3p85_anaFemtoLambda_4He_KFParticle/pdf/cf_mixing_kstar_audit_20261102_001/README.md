# High-k* baseline and mixing/k* audit

- [Figure report (8 pages)](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_mixing_kstar_audit_20261102_001/lambda_he4_mixing_kstar_audit_20261102.pdf)
- [Japanese findings and remaining work](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/Todo/2026/11/02/cf_mixing_kstar_audit_20261102.md)
- [Validation](validation.json)
- [ROOT kinematics test](kstar_prf_math_audit_20261102_result.json)
- [Production function comparison](runtime_functions_verified.json)
- [Exact production diff](runtime_StFemtoMaker_diff.patch)
- [Mixing QA summary](mixing_qa_summary.json)
- [Baseline settings](baseline_settings.yaml)
- [Plot settings](plot_settings.yaml)
- [Numerical audit](baseline_audit.json)
- [Raw/L CF comparison](plots/baseline_group_L.pdf)
- [Fine centralities](plots/baseline_fine_L.pdf)
- [Centrality-composition diagnostic](plots/baseline_cent_exposure.pdf)
- [Species/mass controls](plots/baseline_species_control.pdf)

No physical normalization plateau is established. Kinematics and elementary ME bookkeeping tests pass; the broad baseline origin remains unresolved. Existing ROOT/PDFs and production settings are untouched. `baseline_audit.json` retains all50 CF bins through1 GeV/c; fine fits are reused, not rerun. Species controls use the exact Maker mass gate; L study raw uses the whole-mass-bin gate.

Reproduce the numerical analysis with `audit_mixing_baseline.py --settings <baseline_settings.yaml>`. ROOT probe/drawing commands are recorded in `probe_command.json`, `probe_controls_command.json`, `draw_command.json`, `draw_control_command.json`; they use the mainconf and batch-matched ROOT5 wrapper. Plot CropBox processing preserves numerical data; `plots/crop_manifest.json` records source/output hashes.
