# Λ–⁴He SE/ME k* overlay

This study overlays the pair spectra before CF division. SE remains at its original count scale. ME is multiplied by one constant per centrality group so the SE and scaled-ME integrals agree over **0.4 ≤ k* < 0.5 GeV/c**. Matching the interval integral does not require equality in every bin.

## Input and selection

- 2021 P24iy only; SUMS `1CE33994716F7896423A80AB994BAB14`, verified merge of 26,945/26,994 jobs (49 excluded).
- Mainconf: `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`.
- Config snapshot: `job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14.txt`.
- Selected inclusive ⁴He with the Λ signal mass gate, 1.10708–1.12484 GeV/c². No additional rigidity, nσ, mT or pair selection.
- Original ROOT histograms: `hKstarSEVsCent_lambda_he4_signal` and `hKstarMEVsCent_lambda_he4_signal`.
- Centrality 0–20% = cent9 6–8; 20–80% = cent9 0–5; 0–80% = cent9 0–8. Both SE and ME use the same cent9 range.
- No background subtraction or purity correction. The production Maker's stored-CF normalization defaults are not used here; this explicit diagnostic setting matches the previous high-k study.
- The track-merging QA ROOT is a separate sample and is not an input to this plot.

## Drawing

ROOT 5.34/38 in the batch-matched STAR environment projects each TH2 over the requested centrality range. It sums four native 0.005 GeV/c bins into each displayed 0.020 GeV/c bin, with no smoothing. Only regular k* bins enter the normalization. Each panel establishes the SE axes, shades the normalization interval, then calls `ME->Draw("HIST SAME")` and `SE->Draw("E1 SAME")`. No CF histogram or ratio panel is made.

Black points are raw SE with projected statistical errors. The red line is scaled ME. The vertical axis is pairs per displayed bin, in SE count units. SE point errors do not include event or mixed-pair reuse correlations. Stored scaled-ME errors treat the fitted normalization factor as fixed.

Each page shows 0–1 GeV/c on the left and a 0.3–1 GeV/c detail on the right. Both horizontal and vertical axes are linear. The gray band marks the normalization interval.

| Centrality | SE in norm band | ME in norm band | ME scale |
|---|---:|---:|---:|
| 0–20% | 122,358 | 4,461,431 | 0.027425729547313406 |
| 20–80% | 18,213 | 666,345 | 0.0273326880219706 |
| 0–80% | 140,571 | 5,127,776 | 0.027413638973309288 |

## Artifacts

- [Main 0–20% PNG](se_me_cent0_20.png) / [PDF](se_me_cent0_20.pdf).
- [20–80% PNG](se_me_cent20_80.png) / [PDF](se_me_cent20_80.pdf).
- [0–80% PNG](se_me_cent0_80.png) / [PDF](se_me_cent0_80.pdf).
- [Combined three-page PDF](se_me_overlay.pdf).
- [ROOT histograms and canvases](overlay_histograms.root): native and displayed SE, raw ME and scaled ME for each group, plus the ME scale parameters. Original histograms are preserved.
- [Settings](settings.json), [normalization table](normalization.tsv), [all native/displayed bin contents and errors including flow](spectra.tsv).
- [ROOT validation](validation.json), [independent validation](independent_validation.json), [completion record](completed.json).
- [ROOT commands](root_commands.txt), [ROOT log](root.log), [source](source/DrawSEME.C), [runner](source/run_overlay.py).

The independent reference is the original high-k study's verified `histograms.json`. The physical TH2 projections are compared with all-nσ projections (including nσ underflow/overflow) from the selected-signal QA histograms used by the previous CF. The regular-bin counts must agree exactly. `true/` is a ROOT directory name and does not imply MC truth.

## Reproduction

From the repository root, create a fresh study directory, copy `settings.json` and `source/`, then run the batch-compatible Python with the new settings path:

```sh
/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/python-3.8.13-myixjs7jesskpbhwiowibizezfh52rky/bin/python3 NEW_STUDY/source/run_overlay.py NEW_STUDY/settings.json
```

The runner refuses to overwrite an existing output ROOT. The input ROOT is opened READ, with its recorded size, mtime, inode, device and UUID checked. Fingerprints are checked again after drawing. ROOT output, PNG/PDF signatures, PDF pages, bin-by-bin constant scaling, rebin count/error closure, and normalization closure are verified. No production source/configuration, farm job, mail or Git operation is part of this study.
