# Λ–4He Correlation Function：過去資料と更新ROOTでの再作成準備

> 本書は準備段階の記録。後続の[CF再計算・結果の記録](lambda_he4_cf_update_20261030.md)では試料情報の対応を実装し、計算を開始した。

2026-10-30（日本時間）。対象は2021 P24iy、Au+Au 3.85 GeV固定標的解析。3.85はビームエネルギーのラベルである。

過去の再現用マクロ・README・図付きPDFを特定し、採用済み4He PID/TOFで再解析した統合ROOTの入力検査を完了した。**新CF、Λ purity、LL fitの再計算は未実施**。本計算の前に、レポート生成処理の旧試料固定情報を新試料のmanifestに対応させる必要がある。

## 1. まず参照する過去資料

主要な前回作業は[2026-10-03の日報](../../../../analysisnote/20261003/summary20261003.md)。直接の比較対象も2021 P24iyであり、旧2019解析とは区別する。

| 図付きレポート | 頁数 | 内容 | 説明・再現記録 |
|---|---:|---|---|
| [総合レポート](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/summary_4He_KFParticle_3p85.pdf) | 10 | Λ質量・SE/ME・sideband・背景差引き方式の比較、探索的LL fit | [README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/README.md) |
| [Λ質量窓の低側／高側比較](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves/mass_halves_comparison.pdf) | 13 | full・low・highのCFをRaw/SBsub/CFsub/Method3/Method5で比較 | [README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves/README.md) |
| [mT別CF：0–10%](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10/cf_mt_cent0_10_constant_vs_method5.pdf) | 8 | 分位点4分割＋固定0.2 GeV/c²幅。一定purityとMethod5の比較 | [README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10/README.md) |
| [mT別CF：0–20%](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_20/cf_mt_cent0_20_constant_vs_method5.pdf) | 8 | 同じ比較をcent9=6–8で実施 | [README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_20/README.md) |

PDF 21本の実在、ページ数、サイズ、SHA-256は[資料manifest](assets_manifest.json)へ保存した。総合reportと単独図14本は、既存のreport_metadataに記録されたhash・頁数とも一致した。主要比較reportはそれぞれ上表の13/8/8頁を実測した。

関連する[旧tools/summary_4He_3p85.pdf](../../../../../tools/summary_4He_3p85.pdf)、[tools/summary_4He.pdf](../../../../../tools/summary_4He.pdf)は各10頁。こちらは古い参考資料として残すが、今回その入力来歴を完全に再同定していないため、固定S/Nや設定を新解析へ流用しない。13.5 GeV KFParticle reportは2020 P24iyの別試料であり、今回の対象外。

## 2. 旧入力と今回の入力

| 項目 | 10/3のCF研究 | PID/TOF更新後：今回使う入力 |
|---|---|---|
| 年度・production | 2021 P24iy | 2021 P24iy |
| SUMS ID | `0BEFC64E62F19BBA534A129CABC9B106` | `04BF2DAE47164B6F040E3A460E539289` |
| 採用subjobs | 26,956 / 26,994（38件欠損） | 26,955 / 26,994（39件除外、承認済み部分統合） |
| 入力 | [旧2021 partial_merge ROOT](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root) | [今回の通常ヒスト統合ROOT](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root) |
| PDFの扱い | 旧結果として保持 | 新しい独立出力先へ作成予定 |

今回の正式なディレクトリ名は `auau3p85_anaFemtoLambda_4He_KFParticle`。呼称の`auau3p85GeV_...`と区別する。

完成・統合の根拠は[04BF出力README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/README_outputs_04BF2DAE47164B6F040E3A460E539289.md)、[出力manifest](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_outputs.json)、[投入時runmeta](../../../../../job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289.json)。通常ROOTは14,064,010,756 bytesで、10/26に通常ROOT＋Sparse 3分割の保存・検証・元subjob整理が完了している。

### 今回の実ROOT検査

[検査結果JSON](input_histogram_check.json)と[読み取り専用マクロ](check_cf_input_20261030.C)を保存した。実行環境はSL24y、ROOT 5.34.38、sl73_x8664_gcc485。fitやCF計算は行わず、必要オブジェクトの存在・型・軸・Entriesを読んだ。

- 通常ヒスト数 **1,786**、選択した必要ヒスト **175 / 175** の存在・軸・有限かつ正のEntriesが一致。
- `hRefMult` Entries：**2,100,726,543**、`hRefMultCorr` Entries：**2,067,647,474**。出力manifestの入力／再構成event数と一致。
- `hLambda_InvMass` Entries：**40,118,234**。これはΛ候補のEntriesであり、event数とは区別する。
- `FemtoLambdaRunStatus=completed`、採用4He較正とPID領域のmetadata、Sparse分割manifestを確認。
- 通常ROOTには巨大な`hNuclearPidPairJoint`を含まず、分割manifestを保持。
- 検査対象は全ヒストの読み出し検証ではなく、CFに必要な175ヒストと通常ヒストkey数の確認。各binのCF値・fit成立性・低k*統計の評価は次段階。

初回検査は存在しない1D low/high形式の36キーを仮定してFAILとなった。YAMLと`ExtractMassHalfCF.C`を照合し、実際の保存形式へ検査を訂正した。初回結果は[input_histogram_check_initial_assumed_keys.json](input_histogram_check_initial_assumed_keys.json)に残した。解析ROOT・ヒスト定義を変更していない。

### 必要ヒストの対応

| 用途 | 実際の保存名・数 | 軸・型 |
|---|---|---|
| Λ質量・event QA | `hLambda_InvMass`, `hRefMult`, `hRefMultCorr`：3 | TH1F |
| k*別Λ質量fit | `true/hKstarMass_4He_CentBin0..8`、`mix/hKstarMass_Mixed_4He_CentBin0..8`：18 | TH2F、k*:200bin 0–1、mass:200bin 1.05–1.25 |
| 従来signal/SBのk* | `true/hKstar_4He{,_SBNeg,_SBPos}_CentBin0..8`、mix側は`hKstar_Mixed_4He...`：54 | TH1F、200bin 0–1 |
| full/low/high/SBのk*×cent | `hKstar{SE,ME}VsCent_lambda_he4_{signal,signalLow,signalHigh,leftSB,rightSB}`：10 | TH2F、k*:200bin 0–1、cent9:9bin −0.5–8.5 |
| mT分割 | `true/hKstarMt_4He{,_SBNeg,_SBPos,_signalLow,_signalHigh}_CentBin0..8`、mix側は`hKstarMt_Mixed_4He...`：90 | TH2D、k*:200bin 0–1、mT:200bin 0–10 |

k*・運動量はGeV/c、massとmTはGeV/c²。low/highの1D legacyキーは作られておらず、質量半窓比較は共通のk*×cent TH2を使う。TH2FをTH2Dへ固定castせずTH2として読み込む必要がある。

検査対象は保存Sumw2を持たない。抽出コードは`GetSumw2N()==0`ならunit-fill Poisson（分散=count）、存在すれば保存誤差²を使う。この仮定とTH2F丸めを新reportにも明記する。fitや規格化の相関までこの読み取り検査で保証したものではない。

### 巨大PID Sparse 3分割の位置付け

`..._pid_pair_joint_part000.root`、001、002は入力subjob群10,000 / 10,000 / 6,955件の分割で、RやTOF領域ごとの分割ではない。通常のCFとΛ質量fitは上記の統合ROOTから作成できるため、まずSparseを読む必要はない。

Sparseはk*、原子核R、nσ、cent9、核種表現、SE/ME、排他的Λ質量領域、selected/prePIDの8軸。連続Λ質量やmT軸はない。R別QA射影には有用だが、R別／mT別の独立した連続Λ質量fitをこのSparseから復元することはできない。必要なら各partを射影後に加算し、巨大8Dオブジェクトを再haddしない。β–η等の候補アクセプタンスTH2とペアCFも区別する。

## 3. 採用済み選択条件と投入時設定

入口は[現行mainconf](../../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml)。核PIDは[採用nuclearid YAML](../../../../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_he4_accepted_20261022.yaml)と[4He較正係数](../../../../../config/calibration/dedx_he4_auau3p85_2021P24iy_accepted_20261022.yaml)に記録される。

- 4He：厳密に **0.61 < R < 5 GeV/c**。R<2はTOFを必須とせず、R≥2は有効TOFとdのm²窓 **3.268773–3.693147 (GeV/c²)²** を要求する。境界R=2は高R側。
- d/t/3He：従来Git版のdE/dx関数・係数、既存の種別R下限・vetoを維持し、TOF切替は同じR=2。これらの純度再推定は今回の範囲に含めない。
- 3Heと3He-as-4Heは背景control。後者は同じ3He trackの質量仮説を4Heへ置き換えた表現であり、独立試料や自動的な汚染差引き量ではない。
- 確定した4HeのPID純度vsRを、そのままΛ–4Heペアのpurity(k*)へ代入しない。今回探したΛ質量背景差引きとは別の補正問題である。

[投入時設定照合](submitted_config_check.json)ではruntime bundleのhash一致と12設定を確認した。10設定はbyte単位で一致。mainconfとanalysis_infoの2件はコメント1行ずつの削除だけで、parse後のYAML内容は一致した。投入時analysis_infoは`production_3p85GeV_fixedTarget_2021`、year2021、P24iy、SL24yを指定している。[投入入力リストの照合](dataset_input_check.json)では、保存済みXMLのcatalog URLとSUMSの全26,994 subjob入力リストを確認。81,082ファイルはすべて2021 P24iyパスで、重複なし、subjobリスト連結とmaster listも一致した。リストの付随event値はcatalog metadataであり、実測Entriesと区別する。実際の有効統計は上記部分統合件数とevent数で記載する。この設定照合は、実行バイナリや全解析コードの同一性検証を意味しない。

Λ質量窓は[Maker設定](../../../../../config/maker/maker_auau3p85_anaFemtoLambda_4He_KFParticle_he3Control.yaml)から読む。中心1.11596、σ=0.00296、signal ±3σ、sidebandは左右3–12σ。

| 窓 | Makerの候補質量cut [GeV/c²] | 旧studyのwhole-bin積分 [GeV/c²] |
|---|---|---|
| full | [1.10708, 1.12484] | [1.107, 1.125) |
| low | [1.10708, 1.11596) | [1.107, 1.116) |
| high | [1.11596, 1.12484] | [1.116, 1.125) |
| leftSB | [1.08044, 1.10708) | [1.080, 1.107) |
| rightSB | (1.12484, 1.15148] | [1.125, 1.151) |

保存mass binは1 MeV/c²幅。候補cutとbin全体の積分には境界差があるため、新入力のinspectへ実際のbin範囲を保存する。ペアの運動学的Λ質量仮説は1.115683固定、質量fitの軸は再構成Λ質量である。

## 4. マクロ・README・再現設定

運用入口は[tools/kfparticle_summary/README.md](../../../../../tools/kfparticle_summary/README.md)。原本へ入力ROOT名だけを書き換えて直接実行する手順は使わず、runnerで検査済みコピーを独立出力先へ生成する。

| 処理 | 実装 |
|---|---|
| 一括検査・描画・summary | [run_summary.py](../../../../../tools/kfparticle_summary/run_summary.py) |
| 旧テンプレートを新入力向けに生成 | [prepare_macros.py](../../../../../tools/kfparticle_summary/prepare_macros.py) |
| Λinclusive質量・proxy・診断 | [InspectKFSummary.C](../../../../../tools/kfparticle_summary/InspectKFSummary.C)（fitも行うため、今回の読取検査には使っていない） |
| 監査付きgaus+pol2 / Method5共通処理 | [KFSummarySupport.h](../../../../../tools/kfparticle_summary/KFSummarySupport.h) |
| 総合report・検証 | [build_report.py](../../../../../tools/kfparticle_summary/build_report.py) / [validate_outputs.py](../../../../../tools/kfparticle_summary/validate_outputs.py) |
| full/low/highの比較 | [compare_mass_halves.py](../../../../../tools/kfparticle_summary/compare_mass_halves.py) / [ExtractMassHalfCF.C](../../../../../tools/kfparticle_summary/ExtractMassHalfCF.C) / [report_mass_halves.py](../../../../../tools/kfparticle_summary/report_mass_halves.py) |
| mT別比較 | [run_mt_cf.py](../../../../../tools/kfparticle_summary/run_mt_cf.py) / [ExtractMtCF.C](../../../../../tools/kfparticle_summary/ExtractMtCF.C) / [append_mt_report.py](../../../../../tools/kfparticle_summary/append_mt_report.py) |

原本は[tools/auau3p85GeV/](../../../../../tools/auau3p85GeV)、前回実行したコピーは[旧pdf/macros/](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/macros)。[旧macro_manifest.json](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/macro_manifest.json)に原本と生成物のhash・変更内容が残る。

今回の[互換性照合](macro_compatibility.json)では13原本・13旧生成マクロのhash一致と、旧manifestの条件を使うadapt文字列検査13件の成功を確認した。これはテンプレート互換性の確認であり、新ROOTに対するコンパイル・fit・CF数値動作の検証ではない。適用マクロの新規生成・実行はしていない。

| 元マクロ名（`tools/auau3p85GeV/`） | 出力PDF | 頁 |
|---|---|---:|
| [Draw_LambdaInvMass.C](../../../../../tools/auau3p85GeV/Draw_LambdaInvMass.C) | [Lambda_InvMass_legacy.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/Lambda_InvMass_legacy.pdf) | 1 |
| [draw_LambdaNuclearId_relamom_sum.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom_sum.C) | [kstar_sum_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_sum_4He.pdf) | 12 |
| [draw_LambdaNuclearId_relamom_sum_CFsub.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom_sum_CFsub.C) | [kstar_sum_CFsub_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_sum_CFsub_4He.pdf) | 12 |
| [draw_LambdaNuclearId_relamom_sum_CFsub_method3.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom_sum_CFsub_method3.C) | [kstar_sum_CFsub_method3_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_sum_CFsub_method3_4He.pdf) | 6 |
| [draw_LambdaNuclearId_relamom_sum_CFsub_method5.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom_sum_CFsub_method5.C) | [kstar_sum_CFsub_method5_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_sum_CFsub_method5_4He.pdf) | 6 |
| [draw_LambdaNuclearId_relamom_CFsub_method3.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom_CFsub_method3.C) | [kstar_cent_CFsub_method3_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_cent_CFsub_method3_4He.pdf) | 3 |
| [draw_LambdaNuclearId_relamom_CFsub_method5.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom_CFsub_method5.C) | [kstar_cent_CFsub_method5_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_cent_CFsub_method5_4He.pdf) | 3 |
| [Fit_LambdaNuclearId_relamom_sum.C](../../../../../tools/auau3p85GeV/Fit_LambdaNuclearId_relamom_sum.C) | [Fit_LL_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/Fit_LL_4He.pdf) | 2 |
| [Fit_LambdaNuclearId_relamom_sum_CFsub.C](../../../../../tools/auau3p85GeV/Fit_LambdaNuclearId_relamom_sum_CFsub.C) | [Fit_LL_CFsub_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/Fit_LL_CFsub_4He.pdf) | 2 |
| [Fit_LambdaNuclearId_relamom_sum_CFsub_method3.C](../../../../../tools/auau3p85GeV/Fit_LambdaNuclearId_relamom_sum_CFsub_method3.C) | [Fit_LL_method3_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/Fit_LL_method3_4He.pdf) | 2 |
| [Fit_LambdaNuclearId_relamom_sum_CFsub_method5.C](../../../../../tools/auau3p85GeV/Fit_LambdaNuclearId_relamom_sum_CFsub_method5.C) | [Fit_LL_method5_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/Fit_LL_method5_4He.pdf) | 2 |
| [draw_LambdaNuclearId_relamom.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom.C) | [kstar_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_4He.pdf) | 6 |
| [draw_LambdaNuclearId_relamom_CFsub.C](../../../../../tools/auau3p85GeV/draw_LambdaNuclearId_relamom_CFsub.C) | [kstar_CFsub_4He.pdf](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/kstar_CFsub_4He.pdf) | 6 |

runnerは`nuclearid/`→`true/`、TH2汎用読取、明示ROOT/出力先、4Heのみ、Maker質量窓、中央度群、固定S/Nの置換、rebin・規格化・Method5独立格子、無効fitの扱いをまとめて整合させる。旧toolsの定数（S/N=2.172や5.217）、旧inclusive値0.9672106547576258、旧Method5純度表は新入力へコピーしない。

## 5. 初回の再作成で揃える条件

既存の6 settingsを[settings/](settings/)へ内容を変えず複製した。元settingsとSHA-256を[次回計画](next_run_plan.yaml)に記録する。これらは最初の旧結果比較に使う後処理条件であり、production cutの変更ではない。

| 条件 | 値・扱い |
|---|---|
| 元k* bin | 5 MeV/c（200bin、0–1 GeV/c） |
| CF出力bin | 20 MeV/c、rebin=4 |
| 図のk*範囲 | 0–0.3 GeV/c |
| 規格化範囲 | **0.3–0.4 GeV/c**、出力5bin。表示範囲外も計算へ残す |
| Method5 purity格子 | 50 MeV/c、purity_rebin=10 |
| Method5 fit窓 | 左右2bin＋中心の最大250 MeV/c。端では短縮 |
| massfit | gaus+pol2、1.08–1.15、詳細はsummary_settings.yaml |
| 中央度2群 | 20–80%：cent9 0–5、0–20%：6–8 |
| 中央度3群 | 0–10%：7–8、10–20%：6、20–60%：2–5 |
| 質量半窓比較の追加群 | 0–80%：0–8 |

中央度はSE/MEのraw countを合算してからCFを規格化する。規格化済みCFを単純平均しない。Maker内のQA既定値norm0.5–1.0、rebin5とは別設定である。

### 比較するΛ背景処理

AをsignalのSE、BをsignalのME、L/Rを左右sideband、αを幅補正、N(X)を0.3–0.4でSE/MEを規格化する操作とする。

| 方式 | 定義・特徴 |
|---|---|
| Raw | signal SE/MEを規格化。Λ背景補正なし |
| SBsub | `N[(A−α(SE_L+SE_R))/(B−α(ME_L+ME_R))]`。nominal full α=1/3、half α=1/6 |
| 一定purity CFsub | `(C_signal−(1−PΛ)C_SB)/PΛ`。signal/SBはそれぞれ独立規格化。PΛは当該入力のinclusive proxy |
| Method3 | 各k* binのΛ質量をSE/ME別にgaus+pol2 fitし、窓内のN−Bを求めて比を規格化 |
| Method5 | 独立50 MeV格子の滑動窓fitからP_SE(k*)/P_ME(k*)を求め、元5 MeV countに重み付け後20 MeVへ合算して規格化 |

Method5では40–60 MeVのCF binが50 MeV purity境界を跨ぐ。各元binに対応するPを掛けてから合算し、20 MeV bin中心のPで一括処理しない。inclusive whole-bin proxyの幅比は旧入力で0.339623だったが、新入力でも実際のmass bin境界から計算し記録する。

過去資料は方式比較であり、一つの方式が最終的に優れていると確定した記録はない。数値収束したgaus+pol2でも大χ²・小fit確率が多いため、fitの診断と失敗点を併記する。ゼロMEや不成立fitはNAとし、0や任意のpurityへ置き換えない。

### mTとLL fitの扱い

旧mTレポートの分位点境界2.4/2.6/2.7/2.8/3.9は旧sampleの結果なので、新入力では分位点を再算出する。固定幅比較は2.4/2.6/2.8/3.0/3.2を使う別設定。

mT別Method5は同一centralityのmT積分P_SE(k*)/P_ME(k*)を共有した参考比較であり、mT別独立fitではない。`--shared-method5-approved`が必要な近似であるため、新mT計算を行う段階でこの扱いを明示する。今回mT計算は行っていない。

旧LL fitは0–0.4 GeV/cで規格化範囲と重なる。固定a0=4.01 fm、reff=1.84 fmと全freeの探索比較。描画上限0.3とfit上限0.4を混同しない。規格化、重なるpure fit窓、観測数と背景fit、イベント／混合ペア再利用による共分散の未反映部分があるため、誤差やχ²から最終的な物理結論を出さない。

## 6. 本計算前に行う更新と実行順序

1. `build_report.py`の`SAMPLES['3p85']`は旧SUMS・26,956件に固定され、入力名suffixを`2021P24iy_partial_merge.root`に限定している。新04BFの出力manifestを明示的に読み、年度・coverage・jobid・ROOT fingerprintを検証し、reportへ記録できるようにする。
2. `report_mass_halves.py`も旧SAMPLESと入力exact filenameに固定されるため、同様に修正する。ファイル名変更やcheck削除だけで通さず、2021 P24iy制約と入力一致検証を維持する。
3. shared toolsの変更と旧結果保持を確認後、下記の新しい独立出力先で再実行する。既定`pdf/`や旧`inspect.json`、旧`--resume`を使わない。run_summaryの出力末尾は`pdf`である必要がある。
4. 新入力からinclusive proxyとSE/ME別massfitを再評価。原本から新macrosを生成し、raw/SBsub/CFsub/Method3/Method5を同じ条件で比較する。
5. 必要ならfull/low/high、続いてmT比較へ進む。mT共有P近似は別途明示する。
6. 新reportの全ページで入力jobid・coverage・条件・図とfit診断を確認し、hash/完成state/入力不変を検証する。

**以下は手順案で、上記1–2の実装後に使用する。今回実行していない。** `--settings`等の引数は現行argparseで確認済み。

```bash
python3 tools/kfparticle_summary/run_summary.py \
  --energy 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root \
  --settings Todo/2026/10/30/settings/summary_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf

python3 tools/kfparticle_summary/compare_mass_halves.py \
  --energy 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root \
  --settings Todo/2026/10/30/settings/mass_half_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/mass_halves

python3 tools/kfparticle_summary/validate_outputs.py \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf
```

今回の到達点は「必要資料の特定、入力ヒスト検証、設定複製、修正点と次回手順の整理」。新CF・purity・図付き新reportは次段階。記録は個人`star_local`内に保存し、共有productionコードと過去ROOT/PDFは変更していない。
