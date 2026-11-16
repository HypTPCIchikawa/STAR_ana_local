# Lambda purity: Poisson L with nonnegative constant refits

Current presentation: [L-only report, sources and validation](l_only_report_20261102_001/README.md). The original 24-page report from this folder is a preceding edition containing comparison overlays; the current canonical PDF excludes those overlays. Numerical fit/CF products in this folder are unchanged.

- [Current L-only summary (20 pages, distinct filename)](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/summary_4He_KFParticle_3p85_PoissonL_only_20261102_v2.pdf)
- [Preserved chi-square report (34 pages)](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_poisson_L_constant_20261102_004/summary_4He_KFParticle_3p85.pdf)
- [Japanese results, examples and limits](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/Todo/2026/11/02/poisson_likelihood_20261102.md)
- [Method3 fit atlas, 125 pages](atlas/atlas_method3.pdf)
- [Method5 fit atlas, 50 pages](atlas/atlas_method5.pdf)
- [L vs chi-square CF comparison](comparison/cf_objective_comparison.pdf)
- [Fine Method3 purity comparison](comparison/purity_fine_method3_comparison.pdf)
- [Fine Method5 purity comparison](comparison/purity_fine_method5_comparison.pdf)

The final result is Poisson L plus the requested conditional background model. A negative fitted quadratic background anywhere in the closed mass range triggers a new Gaussian+constant fit, with both slope terms fixed to zero and b>=0. This is data-dependent model selection and its systematic uncertainty is not included in local fit errors. Fit diagnostics record the primary minimum, selected background, status, covariance, boundary and6/4 free parameters. Failed fits remain missing.

## Records

- `poisson_context.json`: input fingerprint, references, mainconf, commands and compiled source hashes.
- `group_extract/`, `fine_extract/`: explicit flat settings and all1,960 extracted spectrum rows.100 are below the fit gate;1,860 independently checked deviances are retained.
- `empty_bin_toy.json`: ROOT5 check of nine empty bins and one unit-count bin.
- `fixed_cases_comparison.json`: the same nine case identities as the preceding report, including numerical trials.
- `fit_objective_comparison.json`, `analysis_validation.json`: exact count matching and independent deviance diagnostics.
- `normalization_settings.yaml`, `normalization_audit.json`, `methods_comparison_selected.json`: both normalization intervals, selected0.40–0.50 CFs and covariance.
- `cases/`, `atlas/`, `comparison/`: final vector PDFs; corresponding `_raw` directories retain originals before CropBox processing.
- `cases_fit_meta.tsv`, `basis_fit_parameters.tsv`, `atlas_tables/`: plotting rows and stable centered-background parameters, including the constant-fallback flag.
- `metadata_correction.json`: inherited one-fit/no-rule-change caveat corrected without changing fit values, counts or covariance. Asymptotic p-values remain in fit JSON; the plotting column is omitted because those p-values are not calibrated for sparse counts.
- `report_content.json`, generated TeX, `poisson_report_metadata.json`, `validation.json`, `publication.json`: report source, page/hash checks and atomic publication/backup.

## Reproduction

Use a new output label and explicit paths from the star-analyzer root. Existing results are protected from overwrite.

1. `prepare_poisson_lambda_study.py --output-dir <newdir> --group-reference <chiSquareGroupJSON> --fine-reference <chiSquareFineJSON> --mainconf <mainconf> --target-report <currentPDF> --center-polynomial --background-seed-floor 0.1 --maximum-function-calls 50000 --negative-background-constant`.
2. Run the saved `commands` in `poisson_context.json` with the batch-matched ROOT wrapper. The target data are2021 P24iy04BF only.
3. `analyze_poisson_lambda_study.py --output-dir <newdir> --previous-cases <fixedCaseJSON> --case-settings <criteriaYAML> --inspection-max 0.5 --raw-likelihood-dir <trial001> --zero-seed-centered-dir <trial002> --before-fallback-dir <trial003>`.
4. Run `audit_normalization_ranges.py --settings <newdir/normalization_settings.yaml>`.
5. `draw_poisson_lambda_study.py --output-dir <newdir> --plot-settings <plotYAML> --chi-square-comparison <sameNormalizationChiSquareCFJSON>`.
6. Run `prepare_poisson_report_content.py`, `build_poisson_lambda_report.py`, then `validate_poisson_lambda_study.py`, each with `--output-dir <newdir>`. Inspect the PDF.
7. `publish_poisson_lambda_report.py --output-dir <newdir> --backup-dir <new unique backup>` preserves the preceding PDF and checks hashes before replacement.

Trials001–003 retain the raw-coordinate L, centered zero-seed L and centered positive-seed L without the final constant fallback. They are diagnostic history, not the adopted result. Source code, settings and small figures are saved to the personal submodule; ROOT/binaries and artifacts above the existing1MiB limit are excluded. No commit or push is performed.
