# Λ–He4 reduced treeの計画作業

ユーザー指示により、Git保存先の訂正後に`anaFemtoLambda_4He_tree`の計画へ復帰した。旧GitHub repo本体の削除は別残件として保持し、この作業を止める条件にはしない。

計画の原本は[Todoの実装計画](../../Todo/2026/11/16/anaFemtoLambda_4He_tree_plan.md)。現行コード、ouraの保存済みproducer/reader、候補treeの不足項目を読み取り、3系統の独立調査で確認した。

- 第1段はevent/candidate/track情報の保存のみ。SE/MEは保存ROOTを再読する第2段。CF処理は外側のマクロに維持する。
- TVector3の物理momentum、raw primary/global rigidity、canonical Float四元運動量を区別する。娘とHeの共有track bank、Λ topology、PID/TOF、helix情報を計画した。
- d/t-only eventのpool入場、処理順・元job境界・UID、EP valid/binの保存が再現性に必要。
- 対象候補と娘だけでは全track/PID前QAを再計算できない。全QAの保存範囲は未回答のため、容量と再現範囲の比較案として明記した。
- 未較正81runは保留。新treeの容量・所要時間は未測定で、histogram pilotから推定しない。

このsessionで解析実装・実行設定作成・build・PicoDst再読・job投入・commit/pushは行わない。保存したのは計画とTodo/日報のみ。新方式のround-trip/closure/容量測定は今後の実装時の受け入れ条件である。

計画原本は設計・技術の2観点で独立レビューし、必須修正なし。関連文書のlocal link13件が存在すること、共有解析の手元差分が作業前と同一であること、両GitのHEAD/indexを変更していないことを確認した。これは計画文書の確認であり、新しいtree実装の動作検証ではない。

2026-11-17追記：[旧GitHub repo本体の削除](../20261117/summary20261117.md)は完了した。tree解析は引き続き計画段階であり、保存範囲・closure・quotaの検討から再開する。

2026-11-17方針更新：ユーザー回答により、初段ROOTへ必要QAヒストとtreeを同梱する方針に確定した。全既存ヒストをtreeのみから再計算する要件は外し、QAのためだけの全track保存を避ける。初段固定QAと後段条件別QAの区別、重複加算防止、容量実測を[計画原本](../../Todo/2026/11/16/anaFemtoLambda_4He_tree_plan.md)へ反映した。前記のA/B未確定という記録はこの回答で更新された。実装・jobは未着手。
