> Superseded by the [Poisson-L study](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/README.md). This directory remains the chi-square reference.

# CF normalization and purity-fit diagnostics

The current summary PDF is updated to 34 pages. The preserved 10-page report is embedded after a reading guide; the 23-page addendum supersedes the stated normalization comparison. Original Fig.13/14 numbering is retained. Original LL parameters are reference results and are not refitted.

- [Updated summary](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_poisson_L_constant_20261102_004/summary_4He_KFParticle_3p85.pdf)
- [Diagnostic addendum](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/summary_4He_KFParticle_3p85_updated_20261102_addendum.pdf)
- [Preserved original PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_norm_diagnostics_20261102_001/summary_4He_KFParticle_3p85.pdf)
- [CF comparison](figures/cf_methods_3_5_6.pdf)
- [Nine mass-fit/normalized-shape examples](figures/normalized_purity_cases.pdf)
- [Raw ratio diagnostic](figures/raw_ratio_normalization.pdf)
- [Japanese results and interpretation](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/Todo/2026/11/02/cf_normalization_diagnostics_20261102.md)

## Numerical records

- `fine_extract/purity_fit_review.json`: 1,260 saved fits, same04BF ROOT and unchanged settings. `fine_extract/extraction_command.json` and `extraction.log` record the ROOT execution.
- `normalization_settings.yaml`, `normalization_audit.json`, `normalization_sensitivity.tsv`: six explicit normalization intervals and sensitivity, selected0.40–0.50.
- `methods_comparison_selected.json`, `cf_points.tsv`, `purity_points.tsv`: selected-reference CF and full shared-purity/covariance results.
- `selected_purity_cases.json`, `all_purity_case_candidates.tsv`, `cases_fit_meta.tsv`, `cases_fit_bins.tsv`: selection criteria, counts, fit diagnostics and displayed bins.
- `cited_sources.json`: original historical and updated source hashes. After publication the preceding updated-summary hash resolves to the backup recorded in `publication.json`.
- `drawing_raw/`, `figures/`, `quoted_figures_raw/`, `quoted_figures/`: original vector drawings and CropBox outputs. Only margins change.
- `report_content.json`, generated TeX, `diagnostic_report_metadata.json`: report source and figure page/hash records.
- `validation.json`, `publication.json`, `diagnostic_source_code.json`: numerical/PDF verification, atomic update/backup, source hashes.

## Reproduction

Run from the star-analyzer root into a fresh output directory. The existing files are protected against overwrite; use a new run label and matching explicit settings. This record contains no production event loop or new job submission.

1. Recreate fine-cent extraction using `fine_extract/extraction_command.json`, the same mainconf and batch-matched ROOT wrapper.
2. Run `python3 -u tools/kfparticle_summary/audit_normalization_ranges.py --settings <normalization_settings.yaml>`.
3. Run `select_purity_diagnostic_cases.py --fine-input <fineJSON> --grouped-input <10/31JSON> --settings <criteriaYAML> --output-dir <newdir>`.
4. Use the saved `draw_cases_command.json`, `draw_cf_command.json`, `draw_raw_ratio_command.json`; run `crop_method6_figures.py` on the raw vector drawings.
5. Run `prepare_purity_diagnostic_report.py --output-dir <newdir>`, then `build_purity_diagnostic_update.py --output-dir <newdir> --content-json <newdir/report_content.json> --sources <newdir/cited_sources.json>`.
6. Run `validate_purity_diagnostic_update.py` before publishing. Verify the new source settings, case paths and cited source/backup mapping when rerunning after a published update.
7. After visual review, `publish_purity_diagnostic_update.py --output-dir <newdir> --backup-dir <new unique backup>` checks hashes and preserves the exact preceding PDF.

The input ROOT is unchanged. No fit gate is relaxed. Shared Method6 purity errors and normalization covariance are retained; event/mixed-pair covariance and model/transfer systematics are not available. Matching purity does not establish a valid mass model. No Git commit or push is performed.
