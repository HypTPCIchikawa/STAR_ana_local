# Poisson L only: current report

All mass fits use ROOT `SLQNR0`, including `L`. All 15 figure pages and every purity/CF value in the 20-page report derive from the adopted L study. The three non-L comparison files and old purity/toy values are excluded from this report.

- [Current summary under a new name (20 pages)](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/summary_4He_KFParticle_3p85_PoissonL_only_20261102_v2.pdf)
- [Preserved preceding mixed-comparison report (24 pages)](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_L_only_20261102_001/summary_4He_KFParticle_3p85.pdf)
- [Report content](report_content.json)
- [Source and figure hashes](poisson_report_metadata.json)
- [Validation](validation.json)
- [Publication](publication.json)

The numerical fits, stored counts, purity/CF tables and original L-only figures are unchanged. Only the report content and explicit L labels change. The previous editions remain archived outside the current report.

Reproduce from the star-analyzer root with a new output directory:

```bash
python3 tools/kfparticle_summary/build_poisson_l_only_report.py --study-dir <004/pdf> --output-dir <new-report-dir> --target-report <preceding-current-pdf> --report-stem <explicit-stem>
```

This validates provenance and builds without publishing. Inspect all pages, then use `publish_poisson_lambda_report.py` with the new output directory and a unique backup directory. ROOT inputs, production jobs, hadd and Git publication are outside this report update.
