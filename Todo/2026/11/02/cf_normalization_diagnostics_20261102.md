> この文書は先行χ²版の記録。現在のL版は[Poisson再評価](poisson_likelihood_20261102.md)。

# Λ–4He CF：規格化域とFig.13・14の質量fit診断

状態：完了。指定の10/30総合PDFを34頁へ更新し、元の10頁PDFを退避保存した。

## 成果物

- [更新済み総合レポート（34頁）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_poisson_L_constant_20261102_004/summary_4He_KFParticle_3p85.pdf)
- [今回の追補だけ（23頁）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/summary_4He_KFParticle_3p85_updated_20261102_addendum.pdf)
- [共通規格化0.40–0.50でのMethod3/5/6 CF比較](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/figures/cf_methods_3_5_6.pdf)
- [9例のSE/ME質量分布・fit・pull・DrawNormalized重ね描き](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/figures/normalized_purity_cases.pdf)
- [規格化域におけるraw SE/ME比](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/figures/raw_ratio_normalization.pdf)
- [元の総合PDF（10頁）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_norm_diagnostics_20261102_001/summary_4He_KFParticle_3p85.pdf)

## 規格化不成立の原因と対応

20–80% Method5のME純度cell0.30–0.35 GeV/cは、k*=0.20–0.45の滑動窓でfitしている。fit mass範囲1.08–1.15 GeV/c²には**1,130,197カウント**、保存mass軸全域には2,034,926カウントがある。status4、isValid=false、covariance status2、χ²/NDF=35,680.414/64≈557.506。低統計とは説明できず、保存fitのピーク/tail付近に大きな構造的pullが見える。Gaussian中心・幅など未収束値は診断表示だけに使用し、純度には採用しない。

この無効な50 MeV純度cellは20 MeV CFの0.30–0.32、0.32–0.34、0.34–0.36に寄与する。旧規格化0.30–0.40がこれらを含むため、全規格化bin必須の判定でCF全体がNAになっていた。

全5統合中央度・Method3/5/6に共通の参照規格化を**0.40–0.50 GeV/cへ変更**し、15件すべての規格化が成立した。同幅100 MeVで、失敗cellと混在CF binを含まない隣接区間である。0.30–0.36の欠測や質量モデル不適合そのものは修復していない。

20–80%のraw mass-gate SE/ME比は0.40–0.50でconstant χ²/NDF=3.155/4、slope=−0.002608±0.007277。ただし0–20%では10.017/4、−0.007428±0.002817であり、全中央度の物理的plateauを確認したという意味ではない。高k*へ範囲を広げると比が低下するため、規格化の参照規約と平坦性の認定を区別した。

| 規格化域 | 20–80% M5のCF scale変化 | 20–80% M6 |
|---|---:|---:|
| 0.40–0.60 | +2.26% | +1.98% |
| 0.50–0.60 | +3.92% | +3.43% |
| 0.50–0.70 | +6.08% | +5.28% |

0.40–0.50を基準とした範囲感度であり、統計誤差には合算していない。6候補の全数値は[規格化比較JSON](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/normalization_audit.json)と[表](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/normalization_sensitivity.tsv)に保存。

## Fig.13・14から選んだ例

9つのfine centralityを同じ04BF ROOTから追加抽出（1,260 fits）。元のFig.13/14は保存して引用し、対応するSE/ME質量分布を赤/青でDrawNormalizedした。Method3では20 MeV cell、Method5では50 MeVの中心cellと最大250 MeVのfit窓を区別した。全9例の1例目は上記の規格化失敗であり、残りは図中の不一致点・一致点・SE fit失敗を選んだもの。

| 例 | Method・中央度・k* cell [GeV/c] | SE/ME fit範囲内count | 主な診断 |
|---|---|---:|---|
| 2 | M3 50–60%, 0.26–0.28 | 27 / 1,489 | P_SE=0.010±0.423、境界。低統計で不安定 |
| 3 | M3 40–50%, 0.14–0.16 | 28 / 1,550 | P_SE=0.185±0.327、低統計で不安定 |
| 4 | M3 0–5%, 0.06–0.08 | 524 / 17,531 | 純度一致、SE χ²/NDF≈0.705。MEにはモデル上の懸念 |
| 5 | M3 50–60%, 0.28–0.30 | 26 / 1,539 | SE fitを棄却 |
| 6 | M5 60–70%, 0.25–0.30 | 42 / 4,246 | P_SE=0.251±0.247、境界。低統計 |
| 7 | M5 50–60%, 0.05–0.10 | 52 / 2,313 | P_SE=0.457±0.159、境界。低統計 |
| 8 | M5 0–5%, 0.25–0.30 | 89,303 / 3,307,029 | 純度・規格化形状は近いがfit χ²/NDF=18.44/599.79 |
| 9 | M5 40–50%, 0.20–0.25 | 1,083 / 48,935 | 純度差0.06315、形状も異なる。ME χ²/NDF≈7.83 |

低統計で赤青の点が離れるという見方を具体例2/3/6/7が支持する。一方、例8は**純度が一致してもfit品質が良いとは限らない**ことを示す。例9や規格化失敗は数十countのSEだけでは説明できず、形状差・mass model・滑動窓内の信号/BG比などの検証が残る。原因は一意に断定していない。Method5の一致例は明示criteriaを満たす別中央度が1件だけだったため、閾値を緩めず不足を記録した。

## データ・誤差・検証

- 元pdf直下の10/3図は旧0BEF sample、更新対象10/30図は04BF sample。両方を区別して引用し、今回の質量診断は04BFに限定した。
- Method6は同じMethod5 ME purity/errorをSE/MEへ共有し、共通推定量の相関を保持する。
- native5 MeVへ重み付け後20 MeVへ合算。fit model/gate、PID、Maker、mixing、元ROOTは変更しない。
- mass model/ME→SE移行系統、未保存のevent/混合pair相関、背景fit–count共分散などは誤差に含まない。normalized-shape χ²は固定面積の記述指標でformal p-valueは付けない。
- 12 tests PASS、736有効CF点と共分散、15規格化、fine10–20%と旧統合10–20%の140 mass rowsの一致、入力ROOT不変、18 PDF productsの頁/hashを確認。34頁を画像化して目視し、LaTeXのはみ出し警告なし。
- 元10頁のCF/LL-fit図は参照。新規格化のCFは追補を優先し、LL fitは再実行していない。

[検証記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/validation.json)／[PDF反映・元版退避記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/publication.json)／[設定・再現手順](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_norm_diagnostics_20261102_001/pdf/README.md)。

追記：高k*の低下と規格化の解釈は[ME/k*監査](cf_mixing_kstar_audit_20261102.md)を参照。0.4–0.5は比較用スケールであり、物理的plateauとして認定していない。保存された数値・PDFは変更しない。
