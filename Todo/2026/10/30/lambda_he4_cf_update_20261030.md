# 更新PID/TOF ROOTによるΛ–4He CFの再作成

2026-10-30（日本時間）。ユーザーの次段階実行指示により、準備記録に続いてCFの再計算を開始した。**計算・図付きreport・成果物検証を完了**。本書を結果・検証の記録として更新する。


## 完成した図付きPDF

- [総合CFレポート（10頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/summary_4He_KFParticle_3p85.pdf)／[README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/README.md)
- [Λ質量窓 full/low/high 比較（13頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/mass_halves/mass_halves_comparison.pdf)／[README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/mass_halves/README.md)
- [新PDF21本の頁数・hashと完了記録](cf_outputs_manifest_20261030.json)／[総合検証](cf_final_validation_20261030.json)／[fit監査集計](cf_fit_audit_summary_20261030.json)

## 入力と実装

2021 P24iy、SUMS `04BF2DAE47164B6F040E3A460E539289`。26,955 / 26,994 subjobsの承認済み部分統合。入力event **2,100,726,543**、再構成event **2,067,647,474**。元ROOTと旧PDFを保持して独立出力先へ生成している。

- [出力ディレクトリ](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf)
- [run_state.json](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/run_state.json)
- [今回のinclusive Λ質量検査](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/inspect.json)
- [検証済みsample情報](cf_sample_provenance.json)
- [準備資料・旧マクロ／README／PDF一覧](lambda_he4_cf_preparation_20261030.md)

`sample_provenance.py`を追加し、明示output manifestと投入runmetaを照合。ROOT実path、SUMS、採用／除外件数、分割入力件数を検証し、hash確認済みruntime bundle内のanalysis_infoから2021 P24iyを確認する。summaryと質量窓reportは新sampleを反映し、旧試料の26,956件やSUMS情報を使わない。

入力／再構成event数は`hRefMult`／`hRefMultCorr` Entriesを優先して検証し、保存TParameter値は別記する。今回両者は一致した。これはCF後処理の変更であり、PID係数・選択条件・Maker・mixing・元解析ROOTを変更していない。

## 計算条件

- CF出力20 MeV/c、表示0–0.3 GeV/c、規格化0.3–0.4 GeV/c。
- Method5純度格子50 MeV/c、左右2binを含む最大250 MeV/c fit窓。P_SE/P_MEは別評価し、元5 MeV/c countsへ重みを掛けてから20 MeV/cへ合算する。
- Λ信号中心1.11596、σ0.00296、signal ±3σ、左右3–12σ sideband。Raw、count-level SBsub、一定purity CFsub、Method3、Method5を比較する。
- 中央度集約は前回と同じ20–80%/0–20%、0–10%/10–20%/20–60%。質量半窓studyでは0–80%も追加。
- full/low/high質量窓比較まで完了。mT共有P近似を用いる計算は今回は行わない。

ここでのgaus+pol2は**Λの不変質量背景を評価する関数**であり、確定した4He dE/dx PID関数・Purity評価を再変更するものではない。

## 更新入力からのΛ inclusive proxy

左右sidebandを実際のbin幅比で補間する推定：

- P = **0.9672107405 ± 0.0000214126**（統計誤差のみ）。
- signal質量窓：**26,632,765**候補、左SB **535,171**、右SB **2,036,120**。
- 実bin幅比 **0.3396226415**、背景見積もり **873268.642**候補。

これはinclusive候補のproxyで、ペア・k*・centrality依存の純度ではない。旧値は流用せず、新ROOTから再評価した。Λ選択が同じため旧inclusive値に近いことと、4He PID変更によるペアCFの変化を区別する。

診断用inclusive gaus+pol2は数値収束したがχ²/NDF = **539308.38 / 64**で、モデル適合は悪い。fit由来のinclusive純度はCF補正へ使わず、sideband proxyを一定purity比較に使う。SE/ME別のMethod3/5もfit診断と無効点を併記する。

## 検証と進捗

- sample provenanceテスト13件、既存13マクロの契約テスト3件、質量窓関連21件がPASS。
- 今回のROOT上でinclusive評価と原本からの新マクロ生成が完了。13マクロの描画・fitを完了。
- 出力hash、入力不変、規格化欠測の全CFへの伝播と該当LL fitの非実行を検証。過去PDF21本のhash不変も確認した。

## 目視検査で発見した規格化処理の修正

初回Method5の20–80%が約0.5へ不自然に下がっていた。MEの0.30–0.35 GeV/c純度fitが不成立となり、一部のnorm binを欠いたまま残りbinの積分で規格化していたことが原因。

`prepare_macros.py`の生成処理に、規格化域の全binで両SE/ME信号推定が必要というguardを追加した。normの一部に欠測がある場合はCF全体をNAとし、図中に理由を表示、LL fitは行わない。Method5は元binの一部pure fitが無効なら出力bin自体も無効にする既存処理を保持している。

同じ生成コード片を使ったROOT5の挙動テストで、正常例、規格化積分は正でも1bin欠損する例、norm外の欠損、負値、非有限誤差を確認しPASSした。生成契約テスト3件も再確認した。影響するMethod3/5の6マクロのみ再生成・再実行を完了し、初回summaryと影響図は`revisions/before_complete_normalization_guard_20261030/`へ保存した。

質量窓比較側には既に全norm bin必須の判定があり、同じ条件ではNAになる。今回の修正でsummary側の扱いも統一する。fitが通るよう閾値や関数を変える操作は行っていない。

## 結果の読み方と残る制約

- Raw/SBsub/一定purity CFsub、Method3/5、中央度別CFと探索的LL fitを更新ROOTから作成した。旧purityや旧CFは流用していない。
- Method5の20–80%は規格化域ME fitが不成立のため**CF全体をNA**とした。9中央度表示でもMethod3のcent9=0,1（70–80%,60–70%）、Method5のcent9=0（70–80%）を同じ理由でNAとした。
- Method3は低k*の一部でmassfit背景が非物理的等となり欠測する。質量窓比較の0≤k*<0.2の積分は必要な全binが揃わず、fullは全中央度群でNA。一方、一部のCF点自体は有効で図に残る。
- 方法の優劣やLLパラメータを最終物理結果とは認定しない。規格化、fit、共有ペア等の未反映の共分散、inclusive proxyの近似、gaus+pol2の大きなχ²はPDF/READMEに記録している。
- 4Heの採用PID・Purity評価は変更していない。今回はΛ質量背景の方法比較で、4He purity vs Rによるペア汚染補正を新たに導入していない。

### 0–80%：0≤k*<0.2 GeV/cの積分CF

正規化は0.3–0.4。下表は信号窓fullのcountを区間積分して規格化した値で、CF点の単純平均ではない。誤差は質量窓比較コードの近似伝播。

| 方式 | 積分CF |
|---|---:|
| Raw | 1.031113 ± 0.008090 |
| SBsub | 1.027649 ± 0.008638 |
| 一定purity CFsub | 1.027943 ± 0.008395 |
| Method3 | NA |
| Method5 | 1.029612 ± 0.008338 |

## 保存

評価コードと小さなPDF/図は既存方針に従って`star_local/snapshots/`へ選択コピーし、他のsnapshot記録を保持した。ROOT、大型数値出力、build成果物はコピーしていない。[保存一覧](cf_private_snapshot_20261030.json)。共有Gitへのcommit/pushは実施していない。

完了時刻：2026-10-30T15:26:34.401861+09:00。
