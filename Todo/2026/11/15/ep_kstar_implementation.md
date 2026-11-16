# EP8各区分k*の保存

**ユーザー指示で中断（2026-11-15 21:01 JST）。承認モードの変更があっても、この実装を自動再開しない。** 専用設定とテストdraftのみ保存済み。production C++本体は未変更、build・新規jobは未実施。

- [x] userのコード修正依頼・仕様を確認。実mixing EPclassのSE/ME、inclusive/signal/L/Rの8TH3D。
- [x] EP8専用mainconf/hist/analysis_infoとEP1対照fixtureを新規準備。既存設定不変の比較とSHAを保存。
- [x] 新規回帰2本をdraft保存（tests/test_femto_epd_kstar_splits.C、runner）。本体未実装に依存し、未コンパイル・未実行。
- [ ] C++でEPclass受け渡し・新8histのfill・型/軸/全本数のinit検証・EP境界metadataを実装。ソース読取の自動承認レビューtimeoutが再試行でも解消せず、手動承認切替を依頼して停止中。
- [ ] batch整合buildとEP8/EP1回帰。全EP合算closure、両方向ME、invalid class、mass window fanout、flow、旧出力不変を確認。
- [ ] 小規模Pico検証・記録更新。新規本番jobは今回未投入。

20:29 JST。[詳細・新設定](../../../../analysisnote/20261115/ep_kstar_implementation.md)。未較正81runは保留。

21:01 JST追記：上記の承認待ち状態から、ユーザーの明示指示による中断へ移行。未完タスクは履歴として残す。再開はユーザーからの新たな指示を待つ。

2026-11-16追記：ユーザーが11/13本番機能維持とdraft退避を選択。準備した設定・テスト・fixture計9ファイルを[私用drafts](../../../../drafts/ep_kstar_20261116_001/README.md)へ移動し、active設定から除外。本体実装は未実施のまま。公開作業は[11/16 Todo](../16/publish_baseline.md)を参照。
