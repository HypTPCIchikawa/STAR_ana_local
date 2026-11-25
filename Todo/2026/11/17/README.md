# 2026-11-17 作業記録

- [旧star_local_submodule repo本体の削除](git_repository_deletion.md)：削除と404確認、STAR_ana_local/public main/ローカル解析の保全、一時認証消去を確認済み。
- [tree解析の実装計画](../16/anaFemtoLambda_4He_tree_plan.md)：初段ROOTへQAヒスト＋treeを同梱し、後段でSE/MEを再計算する。[実装・検証結果](anaFemtoLambda_4He_tree_implementation.md)：108QA＋1186pair closure、EP1/pairY5再解析がPASS。全統計treeは外挿433.35 GiB。旧2AF削除により基礎容量は確保済みで、代表性・本番分割の評価を残す。job・commit/pushは未実施。
- [容量使用量・削除候補調査](disk_capacity_audit.md)：削除前は旧2AFが約2.073 TB（約88%）。ユーザー指示で旧2AFのみ削除済み。ログ/job生成物など他の整理候補は未削除。

- [旧2AF削除とΔφ*–Δη再確認](old2AF_retirement_and_merging_recheck.md)：13,027本を削除済み。比較元PDF/QA ROOTを保持し、新treeから同ヒストを作って傾向を確認するTodoを追加（再確認は未実施）。
- [新tree解析の本番投入](tree_production_20261117.md)：ユーザー指示で較正388run/73,054filesをcompact初段へ。投入前build・fresh closure・XML監査はPASS。SUMS `9037A6F05FAD08FC8E9BD9286C8A38A5`、実4,058jobsの単回tmux投入完了。全入力coverage・記録SHA・完了通知process確認PASS、farmで解析実行中。
