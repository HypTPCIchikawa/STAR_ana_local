# .csh・成功ログの容量整理（2026-11-11）

**最新：整理・独立検査完了。35.399 GB解放、quota87.21%、空き330.257 GB（未確定分控除後）。**

ユーザーが7.7 GeVの追加作業を保留し、手動承認へ切替えたうえで、11/09に未実施だった約35 GBの削除を依頼した。7.7 GeVの新規解析・hadd・再投入は行わない。既存controllerは11/10 11:57 JSTにcompleted_with_issuesで終了済み、今回停止操作は不要。

## 対象

| 対象 | 本数 | 割当GB |
|---|---:|---:|
| A76/1CEの圧縮保存済み生成.csh | 53988 | 14.153 |
| 同入力.list | 53990 | 14.184 |
| 1CEのhadd採用・検証済み成功stdout | 26941 | 7.062 |
| 合計 | 134919 | 35.399 |

対象IDはA76FAF721CDAB36920D53B9E8EEA9230、1CE33994716F7896423A80AB994BAB14。最終ROOT2本、runmeta、設定・joblist・SUMS圧縮記録、1CE除外49 stdoutと代表4 stdout、A76のstdout、未統合2AFの全記録、7.7 GeV全ファイルは保持する。09:24 JSTの独立scheduler照会でqueue0を確認した。

## 実施状況

- [x] .csh/.listの全basename・sizeをtar目録と照合、現tar SHA256がrunmeta検証値と一致。
- [x] 1CE merged26945入力の採用indexへ成功ログを限定し、除外・代表ログを保護。
- [x] 所有者・regular file・許可suffix/dir・指紋を確認する明示削除計画を保存。
- [x] tmux `star_cleanup_logs_20261111_002`、worker PID2472563の実行を確認。
- [x] 全対象を事前・unlink直前に指紋照合し、全134,919ファイルの削除を完了した。
- [x] 予定全ファイル消失、保持ファイル不変、quota低下と独立検査PASSを確認した。
- [x] 完了記録・日報へ実績を追記した。

初回tmux起動（09:26:32 JST）は削除前に終了し、worker・実行marker・logがないことを確認した。修正版ではbashとPythonのパスを明示し、既存tmux内のworker PID2472563とlog作成を確認した。logが見えるまでの時間差で起動直後の確認が失敗したが、既存workerの生存を検出したため重複起動は行っていない。全件指紋照合は通過し、明示一覧に沿って削除中。Wi-Fi切断中もtmux内で処理が継続していたことを再開時に確認した。これは起動時点の記録で、完了結果は末尾参照。

実際の起動は同directoryの `launch_v2.json`、logは `deletion_run_v2.log` を参照。初回の `launch.json` は未実行の履歴として残す。自動再実行は行わない。

[削除計画・本数](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/plan_summary.json)／[実行launch](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/launch.json)／[実行log](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/deletion_run.log)／[独立保護検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/protection_before.json)。圧縮plan SHA256 `5ebeb81a6789d941a7d56bf7055eebf0937e143173b18c52477dcbdf3cfda934`、helper SHA256 `294f47f50725943dc3bd8b00d0e85172d10a5c147fee9ae150b757e6ee462c37`。


## 整理完了 — 2026-11-11T09:55:26.419533+09:00

134,919ファイルの削除と完了検査が終了した。割当ブロック解放量は **35.399139328 GB**、内容の論理サイズは2.560612604 GB。GB=10^9 bytes。GPFSの小ファイルの割当量が内容サイズより大きいため差がある。

| 分類 | 削除数 | 解放割当GB |
|---|---:|---:|
| 圧縮保存済み生成.csh | 53,988 | 14.153 |
| 圧縮保存済み入力.list | 53,990 | 14.184 |
| 1CE統合採用・検証済み成功stdout | 26,941 | 7.062 |
| 合計 | 134,919 | 35.399 |

削除エラー0、予定ファイルの残存0。実行journal全134,919行の件数・論理bytes・解放ブロック・分類別件数をreadbackし、計画・実行結果と一致した。journalは圧縮archiveからのSHA256読戻し一致を確認してから未圧縮の重複コピーを削除した。Wi-Fi切断中もtmux内で処理が続き、再開時に二重起動は行っていない。

独立検査もPASS。A76/1CEの対象.csh/.list残存0。1CEの除外49＋代表4のstdout全53本は前回inventoryの指紋と一致。A76のstdout26,994本を保持。SUMS tar2本のSHA256、runmeta、公開済み3.85/7.7最終ROOTは不変。設定・source/library・その他IDの未統合資料、7.7 GeV全データを保持した。

quota（2026-11-11T09:55:26.419533+09:00）：今回前 **88.58%** → 後 **87.21%**。使用2.257115 TB / soft2.588045 TB。inDoubt 0.673 GB差引き空き **330.257 GB**。quota使用量の純減は35.343 GB。削除ブロック合計との差には監査資料追加・同時刻の割当処理を含むため、同一値とは扱わない。

11/09 ROOT整理202.368876544 GBと今回35.399139328 GBの合計は **237.768015872 GB**。これは対象ファイルの割当ブロックの合計であり、別日間のquota純減とは区別する。

[summary](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/summary.json)／[削除・quota検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/postcheck.json)／[独立検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/independent_postcheck_resume.json)／[削除計画](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/plan_logs.json.gz)／[全実行journal・証跡の圧縮保存](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261111_001/audit_records.tar.gz)。

7.7 GeVについてはユーザー指示で追加作業を保留と記録した。今回、新しい解析・再投入・hadd・QAは行っていない。git commit/pushも行っていない。
