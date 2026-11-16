# Method3・5・6：Λ純度fit点検とCF比較

- [図付きレポート](method6_cf_purity_report_20261031.pdf)
- [CF比較：中央度2分割／3分割](cf_methods_3_5_6.pdf)
- [SE/ME純度比較：規格化域を含む](purity_se_me.pdf)
- [Method3 全fit atlas](atlas_method3.pdf)／[Method5 全fit atlas](atlas_method5.pdf)
- 入力ROOT: [auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root](../../../auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root)
- 2021 P24iy、SUMS `04BF2DAE47164B6F040E3A460E539289`、26,955 / 26,994 subjobs（除外39件）。
- 計算設定: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/Todo/2026/10/31/settings/method6_settings.yaml`。

## 定義と誤差

- Method3：CF binごとのSE/ME Λ質量fitから別々に信号N−Bを求める。
- Method5：SE/MEそれぞれの粗格子・滑動窓fitの純度を元binへ適用してからCF binへ集約。
- Method6：Method5のME純度とその誤差を、同一推定量としてSE/ME両方へ共通適用。
- 同じ純度を独立な2つの誤差として加算しない。単一純度cell内の比では純度誤差が相殺し、格子跨ぎや規格化では依存が残り得る。
- SE/ME同純度は検証すべき仮定。signal/BGの相関や質量形状が異なれば十分な統計でも純度差が生じる。Method6による一致は構成上のもので、仮定の検証にはならない。
- ME→SEへの移行・モデル系統はME fitの統計誤差に含まれない。規格化や共有純度の共分散は保持するが、fit/count、別fit背景間、event/track/mixed-pair再利用、cross-methodの未取得共分散は未反映。
- 規格化域の全binが必要。不成立fitがあればCF全体をNAとし、残りbinだけで規格化しない。ME規格化fitの失敗はMethod6へ移行しても回復しない。
- Method3/5は前回と同じ定義・入力で比較するが、今回の誤差は規格化と共有入力の共分散を扱う近似。同じerror barの再現や方式差の独立有意度は主張しない。

## fit診断（inspection range内）

| 中央度 | Method | sample | usable/total | boundary | min / median observations | median / max χ²/NDF |
|---|---|---|---:|---:|---:|---:|
| 20-80% | 3 | SE | 18/20 | 2 | 3 / 1539 | 1.0924 / 1.8301 |
| 20-80% | 3 | ME | 20/20 | 1 | 139 / 63839 | 4.9929 / 21.036 |
| 20-80% | 5 | SE | 8/8 | 0 | 1780 / 21783 | 3.5357 / 9.4962 |
| 20-80% | 5 | ME | 7/8 | 0 | 69249 / 8.9767e+05 | 60.174 / 557.51 |
| 0-20% | 3 | SE | 19/20 | 1 | 22 / 8464.5 | 2.3334 / 5.9089 |
| 0-20% | 3 | ME | 20/20 | 0 | 937 / 3.085e+05 | 33.133 / 151.25 |
| 0-20% | 5 | SE | 8/8 | 0 | 11633 / 1.1944e+05 | 16.701 / 50.546 |
| 0-20% | 5 | ME | 8/8 | 0 | 4.1669e+05 / 4.3657e+06 | 504.47 / 1770.3 |
| 0-10% | 3 | SE | 19/20 | 0 | 15 / 6599.5 | 2.261 / 4.876 |
| 0-10% | 3 | ME | 20/20 | 0 | 781 / 2.4216e+05 | 27.734 / 123.13 |
| 0-10% | 5 | SE | 8/8 | 0 | 9224 / 92608 | 13.449 / 40.304 |
| 0-10% | 5 | ME | 8/8 | 0 | 3.3284e+05 / 3.4026e+06 | 421.24 / 1434.7 |
| 10-20% | 3 | SE | 18/20 | 1 | 7 / 1865 | 1.263 / 2.7049 |
| 10-20% | 3 | ME | 20/20 | 1 | 156 / 66330 | 6.8322 / 29.495 |
| 10-20% | 5 | SE | 8/8 | 0 | 2409 / 26828 | 4.5925 / 11.381 |
| 10-20% | 5 | ME | 8/8 | 0 | 83846 / 9.6309e+05 | 88.416 / 341.34 |
| 20-60% | 3 | SE | 18/20 | 2 | 3 / 1535.5 | 1.0924 / 1.8139 |
| 20-60% | 3 | ME | 20/20 | 1 | 138 / 63038 | 5.0118 / 22.214 |
| 20-60% | 5 | SE | 8/8 | 0 | 1769 / 21708 | 4.0849 / 9.539 |
| 20-60% | 5 | ME | 8/8 | 0 | 68594 / 8.8784e+05 | 60.146 / 231.51 |

数値収束とモデル適合の良さは別の診断。Method3の低統計fit、Method5の実fit窓の統計量、両sampleのχ²・pull・規格化質量形状をatlasで確認する。中心cellの少なさだけでMethod5のfit全体が低統計とは断定しない。

## 積分CF

| 中央度 | Method3 | Method5 | Method6 |
|---|---:|---:|---:|
| 20-80% | NA | NA | NA |
| 0-20% | NA | 1.0241 ± 0.00879 | 1.0262 ± 0.00855 |
| 0-10% | NA | 1.0246 ± 0.00984 | 1.0276 ± 0.00959 |
| 10-20% | NA | 1.0446 ± 0.02 | 1.0305 ± 0.0191 |
| 20-60% | NA | 1.1015 ± 0.0266 | 1.0774 ± 0.025 |

## 規格化不成立

- 20-80% / Method 5: ME bin 16: invalid contributing mass fit; ME bin 17: invalid contributing mass fit; ME bin 18: invalid contributing mass fit
- 20-80% / Method 6: SE bin 16: invalid contributing mass fit; ME bin 16: invalid contributing mass fit; SE bin 17: invalid contributing mass fit; ME bin 17: invalid contributing mass fit; SE bin 18: invalid contributing mass fit; ME bin 18: invalid contributing mass fit

## レポートに収録した代表fit

- [20-80% / Method 3 / bin 1](fit_method3_cent20_80_001.pdf)
- [20-80% / Method 3 / bin 6](fit_method3_cent20_80_006.pdf)
- [20-80% / Method 3 / bin 16](fit_method3_cent20_80_016.pdf)
- [0-20% / Method 3 / bin 1](fit_method3_cent0_20_001.pdf)
- [0-20% / Method 3 / bin 6](fit_method3_cent0_20_006.pdf)
- [0-20% / Method 3 / bin 16](fit_method3_cent0_20_016.pdf)
- [20-80% / Method 5 / bin 1](fit_method5_cent20_80_001.pdf)
- [20-80% / Method 5 / bin 3](fit_method5_cent20_80_003.pdf)
- [20-80% / Method 5 / bin 7](fit_method5_cent20_80_007.pdf)
- [0-20% / Method 5 / bin 1](fit_method5_cent0_20_001.pdf)
- [0-20% / Method 5 / bin 3](fit_method5_cent0_20_003.pdf)
- [0-20% / Method 5 / bin 7](fit_method5_cent0_20_007.pdf)

## 再生成

```bash
python3 tools/kfparticle_summary/build_method6_report.py \
  --output-dir /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf \
  --settings /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/Todo/2026/10/31/settings/method6_settings.yaml
```

JSON・PDF・設定・生成codeのhashと頁数はreport_metadata.jsonへ保存。元ROOTとproduction条件は変更しない。
