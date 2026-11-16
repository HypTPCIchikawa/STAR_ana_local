# EP・rapidityと高k* CFのPDF

- [x] hadd確定と入力ROOTを確認。
- [x] EP1/8、pairY10/隣接merge、He4 YLab二分、3centralityのSE/MEを投影。
- [x] 固定[0.4,0.5)規格化、normalization共分散、窓CF・形状比・GLSを計算。
- [x] 108closure・独立数値レビュー、180科学図を作成。
- [x] 日本語22ページPDFを組版し、全ページ描画と図表・結論・日本語font・組版ログを検証。
- [x] 完成PDF・ROOT/CSV・provenance/SHAのリンクを解析ノートに記録。入力ROOTの指紋不変。

[解析ノート](../../../../analysisnote/20261115/ep_rapidity_highk_cf_report.md)。主結論はEP8で改善するが全体は1に戻らず、rapidity依存が大きい。81runは保留。

19:02 JST完了。[日本語22ページPDF](../../../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/cf_ep_rapidity_report_ja.pdf)。追加jobや81run較正は行っていない。norm域依存・Y別SE/ME差の原因・高k*拡張は今後の検討として資料に記録。

## Centrality・pair rapidity詳細を取り込む改訂（2026-11-15）

- [x] 元の22ページPDF・TeX・manifest・組版sourceをarchive/20261115_initial22pへSHA確認付きで保存。
- [x] 最新11centrality/484CF系列と独立検証結果を確認し、centrality・pairY・He4二分の科学レビューを保存。
- [x] centrality summaryとpairY×centrality mapのvector図10本を追加し、日本語本文・結論・QA・来歴を更新。
- [x] 新版33ページを組版し、独立科学レビュー・全ページ描画・日本語font・入力不変・組版警告0を検証して同じPDF pathへ反映。

解析コードupdateの中断は維持。新規Pico読取・本番job無し。

21:44 JST更新完了。[33ページ版と改訂記録](../../../../analysisnote/20261115/ep_rapidity_report_update.md)。旧22ページ版はarchive保持。
