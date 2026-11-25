# EP8各binの高k* CFとtree再現性

[日本語PDF・12頁](../../../share/figure/20261120_lambda_ep8_cf_closure_report/lambda_he4_ep8_cf_and_tree_closure_ja.pdf)

完了済み2021 P24iy・EP8/cutなし・depth100の20期間出力を使用。0–60%、norm[0.4,0.5) GeV/cで各EP binを個別に規格化すると、C[0.8,1.0)は8binすべて0.8990–0.9234。CFは窓積分したraw SE/ME比であり、CF点の単純平均ではない。1GeV/cを超える漸近収束は評価しない。0–5%、5–10%、10–20%でも全8binで同窓の値が1未満。相関を含む有意性は評価していない。

旧直接Pico解析とcompact tree再読の再現性は、同run・10,000入力/9,761再構成、depth20・segment resetの同条件で確認。QA108＋pair1186（非空587）、46,336,786セルのcontent/error/Entries/axes/rawSumw2差ゼロ。He4 signal SE13/ME403、H1/24のため、これは再現性の試験であり高k*物理収束判定用ではない。全統計depth100/期間内継続はME条件を変更した別評価。

## 収録

- p1–2：結論、データ・EP/規格化定義。
- p3–6：EP8各binのCF全域・拡大、scaled SE/ME、数値表。
- p7–8：centrality別heatmap・CF。
- p9–10：旧直接解析/tree再読の重ね描き・ゼロ差分・確認範囲。
- p11：独立cent9のCF。p12：再現情報。

## 検証・再生成

日本語フォント・12頁全render目視・文字欠落/はみ出し0。独立数値12113checks、元20ROOT statと進行中設定111pinsの不変を確認。低統計の60–80%は規格化SE=1/0の場合があり、NA/大きな揺らぎを明示。

```bash
python3 rootfile/lambda_ep8_cf_closure_report_20261120_001/build_report.py --figure-directory share/figure/20261120_lambda_ep8_cf_closure_report
```

集計・作図の再現方法は新studyの `EP8_CF_README.md`、確定proofは `report_validation.json` と `ep8_figure_handoff_manifest.json`。ROOT payloadは保存histのみを再読し、eventの再解析・新規job・元source/config・稼働中watcher変更はしていない。
