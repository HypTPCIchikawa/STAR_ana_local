# EP分割・Centrality別の高k* CF 日本語PDF

依頼：高k*でCFが1へ収束するかを、EP分割数とCentralityごとに図付き日本語資料へまとめる。

- [x] 最新EP1/8/16/32の完了集計、EP8各plane詳細、He4Y4/8関連結果の原本SHAを固定。
- [x] 規格化[0.4,0.5)、同EP-valid母集団・深度100・cutなし、窓積分比・表示誤差の意味を確認。
- [x] EP分割数×六10%区間、EP8各plane×Centrality、SEME形状、独立cent9・He4Y関連図を作成。
- [x] 追加依頼：EP8で厳密なHe4YLab<−0.5／>−0.5を選び、k*0–0.3／0–1.0のCFを0–60%・各10%Centralityで作図。Y<−1tailも含め、exact境界除外とUF/OFを照合、低k表示でも規格化[0.4,0.5)を保持。
- [x] 日本語PDFを組版し、数値表の独立監査、全ページ描画・文字・はみ出し・目視を確認。
- [x] 最終PDFと再生成・検証証跡をprivate側へ保全し、解析日誌と索引へ反映。

出力予定：`share/figure/20261122_lambda_ep_centrality_highk_report/lambda_he4_ep_centrality_highk_report_ja.pdf`。作業原本：`rootfile/lambda_ep_centrality_report_20261122_001/`。過去PDFを上書きせず、完了済みhist/数値から図を作成する。高k平均が1に近いことと平坦形状を区別し、有意度・k*>1 GeV/cの推定は扱わない。

## 手動承認で再開

前回は追加CF抽出のROOT起動が自動承認レビュー初回/一度再試行ともtimeoutで未実行となり、待機理由のTodo/日誌追記も同様に未実行だった。ユーザーの手動切替を受け、抽出log/raw/source_before/command/保全JSON/strictmetricsが全て未作成であることと、実際の抽出processが無いことを確認して再開。主10図と既存部DRAFT13頁は作成済み。追加はEP8の厳密He4YLab<−.5/>−.5（Y<−1tail含む、exact境界除外）×0–.3/0–1のCFで、両表示ともH[.4,.5)、共通rebin4。最終16頁PDF・数値と目視検証を仕上げる。原本W/approval_resume.json。

## 最終16ページPDF完成

EP1/8/16/32×Centrality、EP8各plane、He4Y4/8関連結果と、追加指定EP8の厳密YLab<−0.5/>−0.5 CFをk*0–0.3/0–1.0で掲載。両表示とも規格化H=[0.4,0.5)、rebin4。Y<−1も下側に含め、exact境界/invalidYの残差は全cent・kflow込みで0。

0–60%の高k窓R[0.6,1)はY<−0.5で0.977724、Y>−0.5で0.910332。下側の二窓R[0.6,0.8)→R[0.8,1)は0.983993→0.971853。低Yで1に近づくが平坦さ・有意度を断定せず、EP16/32とのY条件組合せは今回未検証。背景差引き前、alpha固定誤差。

PDF：`share/figure/20261122_lambda_ep_centrality_highk_report/lambda_he4_ep_centrality_highk_report_ja.pdf`。全16頁decode/render、日本語抽出・TeX欠字/はみ出し・数値入力/全図SHA・目視確認PASS。前回のROOT起動と待機記録の保存は承認レビューtimeoutで未実行だったが、手動切替後にpartial無しを確認して抽出を一度実行し完成。本番ROOT/source/lib/config変更・新job投入なし。`delivery_verified.json`とprivate verificationへ証跡保全。
