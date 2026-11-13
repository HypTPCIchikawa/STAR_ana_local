# 2026-11-13 作業記録

**最新更新（2026-11-14 JST）**：[quota整理と個人submodule公開](quota_cleanup_and_local_push.md)。指定7.7 GeV clone・旧テストROOT84本・旧ログ181,376件の削除は完了。圧縮ログと必要な基準データは保持し、個人submoduleへのpushとremote HEAD一致の確認まで完了。

- [Λ–⁴He EP・YLab・track-merging解析コード修正](../../../../analysisnote/20261113/summary20261113.md)：実装とbatch build／回帰PASS、実Pico最終closure確認中。
- [EP較正のrun coverage](../../../../rootfile/femto_lambda_epd_ylab_merging_validation_20261112_001/epd_run_coverage.json)：既存全投入snapshotの約90.2%のfilesに対応。残り81runsの追加較正が必要。

- [コード反映・全検証完了／残較正Todo](../../../../analysisnote/20261113/lambda_he4_epd_ylab_merging_implementation.md)：noEP既存1675hist完全一致、EP1/EP8の1万event比較と3回帰PASS。本番再投入なし。


## EPの全run対応方針の補足 — 10:49 JST / 2026-11-13

[既存係数の共有可否・Pico内情報・EP確認図](../../../../analysisnote/20261113/ep_calibration_transfer_and_all_runs.md)。未対応81runは前期day121–125。今のrun別EP表に係数が無いという状態であり、81runすべて個別に新規較正する以外の方法が無いと断定しない。

- [ ] 対象81runで明示した共有係数候補を独立標本で検証する。
- [ ] 共有が適合しなければ前期の期間共通／run別recenterと、それに対応するshiftを準備・独立検証する。
- [ ] 全runのEP1/EP8再処理。今回source/configの有効化・SUMS投入は未実施。
- [x] 補正後EP角・実8区分占有率を確認する境界線入りPNG/PDFを作成（1万event pilot、有効9374）。

- [81run保留方針・pair rapidity10区分・quota監査](pair_rapidity_and_quota_20261113.md)：実装・検証中。

- [pair rapidity10区分実装・全検証・quota監査完了](../../../../analysisnote/20261113/pair_rapidity_implementation_and_quota.md)：25files/jobの388run用XMLまで準備。本番未投入、81runは保留。

- [22:22 JST 本番投入開始・メール監視の連結](../../../../analysisnote/20261113/pair_rapidity_production.md)：較正済388runだけのEP1/EP8、実SUMS分割で容量再確認。過去の未投入記録よりこちらが最新。81run保留を維持。

- [22:33 JST 本番投入・完了メール監視の起動確認済み](../../../../analysisnote/20261113/pair_rapidity_production.md)：各4,058jobs、全group一致。初期入力読取エラーの回復確認をTodoへ追加。

- [23:01 JST quota削減・star_local push](quota_cleanup_and_local_push.md)：指定7.7 GeV clone削除済み、約25.876GB削減。残る旧ログ/test ROOT整理とpushは承認基盤の手動切替待ち。
