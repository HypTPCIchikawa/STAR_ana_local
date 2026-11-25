# EP8各binの高k* CFと旧解析/tree再現性の資料

ユーザー依頼：進行中EP1/8/16/32 jobの完了を待たず、既存EP8二段階目出力の各binでCFの1への収束を比較し、treeを介した解析と旧解析の一致確認と合わせた日本語の図付きPDFを作る。

- [x] 完了済みEP8 depth100・cutなし20出力を対象とし、raw signal/規格化0.4–0.5、表示上限1 GeV/cを確認。
- [x] 実ROOTを再読し、同10k入力/9761再構成、depth20/reset=trueの旧直接Picoとcompact tree/replayのQA108+pair1186を全bin/content/error/Entries/axis/rawSumw2一致で確認。He4 signal SE13/ME403、規格化域SE1/ME24。
- [x] EP8各binのCF/scaled SE・ME、高k窓積分比、centrality別比較を独立に集計・照合。
- [x] 日本語PDFを作成し、図・日本語フォント・数値/説明の一致を検証。
- [x] final PDFと再現用helper/小容量証跡を個人領域へ保管し、解析記録を追記。

作業原本 `rootfile/lambda_ep8_cf_closure_report_20261120_001/`、資料先 `share/figure/20261120_lambda_ep8_cf_closure_report/`。既存/root/pins/新規本番jobは変更しない。旧直接解析との一致確認の条件と、全統計の期間内継続/depth100へのME条件変更を区別する。

## 完了

日本語PDF12頁を作成・全頁目視確認。0–60%各EP binのC[0.8,1)は0.8990–0.9234、全8binで高k*の1からのずれが残る。0–5/5–10/10–20%でも同傾向。物理起源・漸近値・有意性は断定せず、1GeV上限と規格化/共有event相関、末梢低統計を明記。

[最終PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/share/figure/20261120_lambda_ep8_cf_closure_report/lambda_he4_ep8_cf_and_tree_closure_ja.pdf)

旧新同10k試験は108QA+1186pair完全一致。全統計depth100/期間内継続との条件差を説明し、同一sampleのclosureと全統計物理形状比較を分けた。旧JSON読取のレビューtimeoutは、手動設定再適用後の独立実ROOT監査・新集計で解消。元20ROOT/111pins不変。新規event再処理/submit/source変更0。証跡は `star_local/verification/ep8_cf_closure_report_20261120_001/`。
