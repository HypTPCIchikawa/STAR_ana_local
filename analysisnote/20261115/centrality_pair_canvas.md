# Centrality別pair rapidity SE・ME・CFのROOT直接PDF

## Session update — 21:12 JST / 2026-11-15：完成

ユーザー依頼により、従来の0–20%・pairY10のROOT直接PDFをcentrality別に拡張した。[全体PDF121ページ](../../../share/figure/20261115_lambda_centrality_pair_canvases/centrality_pair_rapidity10_SE_ME_CF_root.pdf)と[centrality別分冊一覧](../../../share/figure/20261115_lambda_centrality_pair_canvases/README.md)を保存。旧10ページPDFは保持している。

対象はcent9全9区分（0–5,5–10,10–20,20–30,30–40,40–50,50–60,60–70,70–80%）と、比較用0–10/0–20%の11slice。各sliceにinclusive＋pairY10の11ページ、計121ページ。分冊11本は各11ページ。ROOT Canvas Divide(2,3)で、左EP1・右EP8、上SE＋scaled ME・中CF full・下CF high。各図にcentrality/cent9/pairY条件、SE/ME総ペア積分、H内ペア数、MEscale、有効CF点数、空bin数を記載した。件数はevent数ではない。

条件は従来を継承：raw Lambda-He4 signal CF、背景差引・purity補正前、norm H=[0.4,0.5) GeV/c、rebin4、k*表示0–1・high0.5–1。pair YLab[-1.5,0.5)の10区分。各EP×cent×Yで別のSH/MHを使い、合算centralityはSE/ME加算後に規格化。共有normのJacobian誤差を継承し、event再利用とEP1/8間相関は未評価。規格化不能・点のSE/ME空binを0/1で埋めない。

## 検証

- ROOT5.34/38のbatch整合wrapperで抽出484系列、SE/ME968ヒストを保存。396投影closureがregular範囲でPASS、最大count差0。
- 0–10=sum0–5＋5–10、0–20=sum0–5＋5–10＋10–20のnative内容・分散（k*flows含む）176比較を独立検証しexact。
- 旧3centralityの264ヒスト・132CF系列（点/誤差/窓/GLS含む）が完全一致。
- 全484系列・24,200点のCF/NA/誤差、484共分散行列、542 GLSを独立再計算し数値精度内で一致。metrics484行も一致。
- 実ROOT描画時の242列のpair数・H count・有効CF点数をJSONと照合し、combinedと分冊の双方で一致。
- 全12 PDFのページ数（121＋11×11）、全121ページ描画、text抽出、全体一覧と代表ページを目視。原merge ROOTのsize/mtimeNs/inode/deviceは開始時から不変。元ROOTの全file SHAは計算していない。
- 9single sliceの全合算はcoverage/totalsを記録したが、別の0–80%独立投影との比較は未実施。peripheral・一部Yの低統計は資料にそのまま残す。

## 再現と成果物

studyは`rootfile/femto_lambda_cent_pair_canvas_20261115_001/`。`study.yaml`、`input_manifest.json`、`projection_plan.tsv`、`extract.log`、`analysis_results.json`、`metrics.csv`、`cf_covariances.npz`、`independent_validation.json`を保存。[CF/raw/scaled hist ROOT](../../../rootfile/femto_lambda_cent_pair_canvas_20261115_001/cf_histograms.root)、[全Canvas ROOT](../../../rootfile/femto_lambda_cent_pair_canvas_20261115_001/canvases_all.root)、[最終manifest](../../../share/figure/20261115_lambda_centrality_pair_canvases/canvas_manifest.json)。作図macroは同studyの`draw_split_canvases.C`と`draw_canvas_batch.C`、実行logは`draw.log`。

既存のmainconfは`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_covered25.yaml`。jobidはEP1 `FD44A5EB1B0AA22C546F169E86482587`、EP8 `DA0DEB4347BD834F1096F8D8798F2133`。configlogは`job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_<jobid>.txt`。2021 P24iy、両系列共通成功4023group・72424 Pico参照。今回の新規mainconf/jobid/configlogは無し。

## 中断方針

[EP別k*保存の解析コードupdate](ep_kstar_implementation.md)はユーザー指示どおり中断。専用設定・未実行テストdraftを保持したまま、production C++は変更していない。今回変更したのは新規作図用postprocessingと記録だけ。新規Pico読取・本番job・メール・pushは無し。未較正81runも保留し、承認モードの変更だけでは解析コードupdateを再開しない。
