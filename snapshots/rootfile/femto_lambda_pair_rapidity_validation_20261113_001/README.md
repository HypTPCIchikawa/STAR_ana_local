# Pair rapidity10区分・81run保留・quota監査の検証記録

[実装仕様と結論](../../analysisnote/20261113/pair_rapidity_implementation_and_quota.md)。

- `build.log` / `build_command.json`: SL24y/gcc4.8.5/ROOT5.34.38のbatch-matched build。
- `pair_unit_result.json` / `pair_unit_retry.log`: 1603559チェックPASS、うちpairY1140649。10区分/flow/SE-ME/全mass/cent/隣接merge/He3除外/既存診断不変/TH3D精度/config拒否/実2workerhadd metadataを確認。
- `pair_projection_unit_result.json`: 実QA helperの130チェックPASS。native ROOTの(0,0)・表示zoomの影響を避け、counts/variance/Xflow/variableedges/source range保持を確認。
- `real_pico_closure.json`: 同じ10000eventsで旧1707physics hist・40002230cells完全一致、Y区分とflowの24240cells閉包。He4信号SE13/ME403不変。
- `pairY10_qa_validated.pdf`: 最終33頁。`qa_validated_pdf_check.json`は全頁Ghostscriptとsignal10区分のSE/ME4頁目視の記録。初期のほかのQA PDFsは中間作業用で、最終判定には使わない。
- `quota_resume_check.json`: 再開後の最新quota。モデルは`quota_projection_models.json`/`quota_review.md`。25files/job・2系列+haddの約166.9GB、ROOT2倍stress219.2GBが最新soft余裕約286.2GB以内。
- `covered388_pico_files.list` / `covered388_input_manifest.json`: 同じ2021P24iy元masterの順序を保持、73054files/388runs、81runs除外。EP1/8で同じinputと実groupingを使う。
- `generate_joblist_commands.json` / `batch_xml_check.json`: 25files/job、約2923jobs/系列のXML生成・検証。**star-submit未実行、本番job未投入**。
- `implementation_manifest.json`: 最終source/config/QA/test/library/XML SHA。
- `changes_from_preexisting_worktree.json`: 編集前14filesのbackupとの比較。既存dirty作業・過去ROOT/snapshotは保持。

未較正81runの解析・較正・共有試験はユーザー指示で保留。EP分割の有用性確認後に較正を検討する。投影・実データ試験はコードの動作確認であり、高k*CF問題の解消を実証した結果ではない。

一時中断後に再開し、最終QA目視・独立source review・差分・manifest・quota再確認を完了した。実投入の際は最新quota、実SUMS job数/file grouping、25fileworkerのruntime/stagingを改めて確認する。実packed25workerの出力全量を測った容量保証ではなく、保守modelを含む予測である。
