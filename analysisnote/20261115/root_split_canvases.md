# EP8各区分・pair rapidity各区分のROOT直接Canvas

## 19:49 JST / 2026-11-15

ユーザー依頼は、EP8各planeとpair rapidity各区分について、scaled SE/MEの重ね描きとCFを条件名付きのROOT Divide Canvasで別PDFにすること。

## Pair rapidity：作成完了

[ROOT直接PDF：全10ページ](../../../share/figure/20261115_lambda_split_canvases/pair_rapidity10_SE_ME_CF_root.pdf)。1区分1ページ、Canvas Divide(2,3)。左EP1・右EP8、上段SE＋scaled ME、中段CF全範囲、下段0.5–1 GeV/c拡大。条件は前回と同じ0–20%（cent9 6–8）、pair YLab[-1.5,0.5)の10区分、norm[0.4,0.5)、rebin4、raw signal CF。

前回の`cf_histograms.root`から、各EP/各Yの`SE`、`ME_scaled`、`CF_graph`を直接cloneして描画。MEの二重scaleやTH1 DivideによるCF誤差の再計算は行わない。共有規格化誤差・欠損点のmaskをそのまま継承。規格化SEが無いbin01/09/10はNA表示、bin02はS_H=10の低統計を明記した。

- [Canvas保存ROOT](../../../rootfile/femto_lambda_split_canvas_20261115_001/pair_rapidity10_canvases.root)
- [全ページ一覧](../../../share/figure/20261115_lambda_split_canvases/pair_contact_sheet.png)
- [設定・入力/出力SHA・実行command](../../../share/figure/20261115_lambda_split_canvases/canvas_manifest.json)
- [描画macro](../../../rootfile/femto_lambda_split_canvas_20261115_001/draw_split_canvases.C)
- [実行log](../../../rootfile/femto_lambda_split_canvas_20261115_001/draw.log)

ROOT5.34/38のbatch整合環境で出力、RETURN_CODE0、DRAW_SPLIT_CANVASES_OK pages=10。Ghostscriptで10ページ確認・全ページ描画。bin02/04を個別に目視、全ページ一覧を確認。保存済み入力ROOT不変。独立レビューで元20系列のraw/CF入力・有効点数とNA/低統計条件を照合した。

## EP8各plane：保存情報不足で未作成

本番EP8 mergeをREAD-only監査。EP8はcent9×Vz×EP8のME pool条件として適用されているが、k*の保存先は全EP合算。`hKstar{SE,ME}VsCent_lambda_he4_signal`はk*×cent9、TH3はpairYの8本のみで、EP joint軸は無い。`hEpdPsi1VsCent`、`hEpdMixingBinVsCent`はevent数のQAであり、k*分布を分割する情報を持たない。

`true`（846 keys）・`mix`（801 keys）のdirectoryにもEP hist/TH3/THn/Treeは無し。`eventLedger`にEP枝がなく、`selectedLambda`はΛ情報だけでHe4情報無し。writeDetailedLedger=false、candidate store disabledで、保存済み候補からの再pairingも不可。従って要求された8区分別SE/ME/CFのPDFは現在の出力から作れない。全体CFを8区分分コピーしたり、event数比で分配したりしない。

必要な追加は、例えば`hKstar{SE,ME}VsEpBinVsCent_lambda_he4_signal`（k*×EPbin×cent9）を実際のSE/ME充填箇所で保存し、その出力を再生成すること。今回の作図依頼ではproduction source変更・新規job・Pico再読を開始していない。この部分は未完了としてTodoに残す。

## 来歴

mainconf：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_covered25.yaml`。jobid EP1=`FD44A5EB1B0AA22C546F169E86482587`、EP8=`DA0DEB4347BD834F1096F8D8798F2133`。configlogは`job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_<jobid>.txt`。詳細入力指紋は[前回記録](ep_rapidity_highk_cf_report.md)と`input_manifest.json`。2021 P24iy、両系列共通4023group、未較正81run保留。今回の追加mainconf/jobid/configlog・メール送信は無し。
