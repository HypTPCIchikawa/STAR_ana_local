# EP分割・rapidity選択による高k* CFの比較

対象は2021 P24iy、Λ–⁴He signal窓のraw CF。haddは2026-11-14 19:03 JSTに完了。共通成功4023 group（72424 Pico references）、未出力35group除外、未較正81runは保留。

## 解析条件

- study: `rootfile/femto_lambda_ep_rapidity_cf_20261115_001/study.yaml`。
- 0–20%（cent9 6–8）を主結果、0–10%/10–20%はcross-check。
- 固定規格化[0.4,0.5) GeV/c、200bin→rebin4で50bin。k*上端1 GeV/c。
- PairY10とHe4 YLab二分はSE/ME共通のpair選択。EPだけがpoolを分割し、rapidity pool変更やjoint二軸選択は今回の出力に無い。
- signal raw CF、背景差引・purity補正前。共有規格化のJacobian covarianceを含め、event/pair再利用とEP1/8相関は未評価。

## 数値結果（PDF反映済み）

0–20% inclusive C[0.8,1)=0.866149±0.003064(EP1)→0.906265±0.003204(EP8)。形状比T=CF[0.8,1)/CF[0.5,0.7)=0.907613→0.935555。EP8で改善するが、1未満と下降が残る。

EP8のHe4 YLab<−0.5は0.968207±0.005376、>−0.5は0.877709±0.004045。pairY[-1.1,-0.7)の隣接mergeは0.981327±0.008475/T=0.994836±0.006068でほぼ平坦だが中心値約2%低い。10–20% He4Belowでは1.00394±0.01560と1・平坦に整合するが全体解決とは区別する。空bin/規格化SE10の低統計binは判定不可。

## 検証

132 CF系列、108投影closureがregular k*範囲で完全一致。独立レビューで全CF・窓積分・誤差・共分散・184 GLSを再計算し一致。k*overflowのMEにEP1=3、EP8=2 countの差があり、対象範囲外として記録。ROOT5で投影264ヒストと180科学図を作成済み、vector crop後に日本語PDFへ組版する。

訂正：初稿のoverflow差についての「legacy float差」という原因の記述は、このhistについて直接検証していないため撤回する。確定している事実は上記のcount差と、CFに使うregular範囲の差が0であること。

## 実行・再現

- 入力ROOT/metadata/mainconfは`input_manifest.json`、投影選択は`projection_plan.tsv`。
- 元mainconf/SUMS ID/configlogは前日の[hadd記録](../20261114/pair_rapidity_hadd_20261114.md)を参照。
- ROOT5抽出：`extract_cf_inputs.C`、実行log `extract_retry.log`。
- 数値：`/usr/bin/python3 -I .../analyze_cf.py`（system NumPy1.23.5）。通常CVMFS Pythonではlibgfortran.so.3が不足するため-isolated system Pythonを使用。
- 描画：`make_cf_plots.C`、log `plots.log`。
- 新規SUMS job/Pico再読/較正/解析source変更/メール送信なし。

---

## Session update — 19:02 JST / 2026-11-15：日本語PDF完成

[比較PDF：日本語22ページ](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/cf_ep_rapidity_report_ja.pdf)を作成した。SEと規格化MEの重ね描き、全体のCF、He4単体Y二分、pairY全10区分・隣接2区分ずつの合算、centrality cross-check、誤差と解釈上の限界を掲載。180本のvector科学図を保存した。

主結論は、EP8で改善するが0–20%全体の高k* CFは未解決。0.8–1.0 GeV/cの窓CFは0.866→0.906で、EP8でも1より約9.4%低く下降が残る。He4 YLab<−0.5では0.968、pairY[-1.1,-0.7)では0.981でほぼ平坦となるが、全体の解決や漸近収束の証明ではない。固定規格化[0.4,0.5)が物理的plateauであるとは未確認。ヒスト上端1 GeV/cのため、評価範囲は0.5–1 GeV/cに限る。

最終検証：PDF22ページを全描画し、一覧と主要ページを目視。日本語text抽出・font埋込、組版のMissing character/Overfull/Underfull/Warning各0を確認。入力ROOTのsize/mtimeNs/inode/deviceは解析開始時から不変。108closureはregular範囲差0。独立数値レビューは132系列と184 GLSに対して実施済み。最終PDF自体の独立agent読み取りは承認レビューtimeoutで未実施のため、その部分はrootによる実行済み検証を根拠とする。

組版時にXeLaTeX擬似boldの描画状態が取り込んだROOT図へ継承される問題を発見し、図の取り込み箇所だけtext rendering modeをresetした。元180図・crop180図・数値・入力ROOTは変更していない。

成果物：

- [derived ROOT：raw SE/ME、scaled ME、CF graph](../../../rootfile/femto_lambda_ep_rapidity_cf_20261115_001/cf_histograms.root)
- [132系列の数値CSV](../../../rootfile/femto_lambda_ep_rapidity_cf_20261115_001/metrics.csv)
- [設定YAML](../../../rootfile/femto_lambda_ep_rapidity_cf_20261115_001/study.yaml)
- [全数値・共分散の参照JSON](../../../rootfile/femto_lambda_ep_rapidity_cf_20261115_001/analysis_results.json)
- [最終provenance・SHA・検証manifest](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/report_manifest.json)
- [全ページ一覧](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/contact_sheet.png)

次の検討は、同じラピディティ条件でのSE/ME形状差の原因、norm領域依存、必要に応じた高k*範囲拡張とevent再利用相関の評価。今回の結果だけでrapidity pool分割が有効と結論しない。未較正81runはユーザー方針どおり保留し、追加較正や新規本番jobは開始しない。

---

## Session update — 21:44 JST / 2026-11-15：中央度別結果を反映した33ページ版

ユーザー依頼により[同じPDF path](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/cf_ep_rapidity_report_ja.pdf)を33ページ版へ更新した。上の22ページ版の作成記録は履歴として保持し、[旧版PDF](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/archive/20261115_initial22p/cf_ep_rapidity_report_ja.pdf)も保存。新しい単独9中央度、pairY×中央度map、He4二分の中央度依存・支持統計を追加。0–30% inclusiveの残留、限定Y条件の整合性、周辺の精度不足を区別して結論を更新した。[改訂記録・検証・来歴](ep_rapidity_report_update.md)。解析コードupdateの中断・81run保留は維持。
