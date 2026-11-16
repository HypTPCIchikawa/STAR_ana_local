# 容量整理（2026-11-09）

ユーザー指示：不要なlog、元作業領域 `job/run/` の生成.csh、検証済みhadd入力・不要なテストROOTを削除する。最終統合ROOTとQA、解析source/config、再現用圧縮記録、実行中7.7 GeV再投入jobが参照するファイルは保持する。

## 作業

- [x] 元作業領域の規約・前回11/04の整理記録を確認。
- [ ] 完了済みjobの.csh/logを、終端状態と圧縮済みSUMS記録に照合する。
- [ ] 3.85 GeV `1CE33994716F7896423A80AB994BAB14` の検証済み統合ROOTを保持し、重複した元ROOT・中間ROOTを一覧化する。
- [ ] 7.7 GeV旧部分統合 `36B7E875A19554DB711364BE0B072DA5` の最終ROOTを保持し、元2090本・中間ROOTを一覧化する。再投入分は保持する。
- [ ] テスト完了証拠のある不要ROOTだけを対象に加える。
- [ ] 削除直前にpath・device/inode・size・mtimeを再確認し、削除一覧と結果を保存する。
- [ ] 削除後に最終ROOT不変・再投入監視の稼働・quotaを確認し、実解放容量をまとめる。

## 保持する未統合データ

`2AF731F2C75544890B4E858354DDE43D` の候補保存ROOT約2.07 TBは、検証済みhaddの重複入力とは確認できていないため削除しない。旧 `auau3p85_anaLambdaNuclearId` の対応未確認ROOTも保持する。

## 状態

削除は未実施。監査担当ごとに限定した候補一覧を `/tmp/yichikawa/cleanup_20261109_001/` に作成中。削除前quotaの再照会・プロセス一覧照会は自動承認review timeoutのため未実行（同一操作の再試行1回でも未実行）。14:44 JSTの実測記録では使用96.10%、inDoubt差引き空き94.45 GB。削除容量はファイルの割当ブロックとhardlinkを区別して記録し、quotaの変動と照合する。


## ROOT整理完了 — 2026-11-09T15:33:48.567318+09:00

検証済みROOT **29,340本**を削除し、最後のhardlinkの割当ブロックで **202.369 GB**を解放した。削除エラー0、予定ファイル残存0。

| 対象 | 件数 | 解放割当GB |
|---|---:|---:|
| 3.85 GeV 1CE統合済み元ROOT | 26945 | 163.716 |
| 同hadd中間ROOT | 273 | 25.729 |
| 7.7 GeV旧部分統合の元ROOT | 2090 | 11.161 |
| 同hadd中間ROOT | 22 | 1.316 |
| 検証済みテストROOT・入力cache | 10 | 0.446 |

GB=10^9 bytes。7.7 stage最終ROOTのhardlink別名削除分は容量解放に数えていない。3.85のstage最終hardlinkは保持。公開済み3.85/7.7最終ROOT2本はsize/mtimeNs/inode/device不変で、ROOT magicとheader fEND=file sizeも確認した。QA PDF・図、source/config/library、全再投入出力を保持した。

削除後quotaは使用 **2.285/2.588 TB（88.30%）**。inDoubt 5.847 GB差引き空き **297.079 GB**。削除前直近の14:44 JSTは96.10%・空き94.45 GBだった。基準時刻が異なり再投入等の書込みもあるので、quota差と202.369 GBの削除ブロック合計を同一とは扱わない。

15:33:35 JSTの7.7再投入controllerはmonitoringで生存。batch10は161終了・39 running、idle/held0。解析やschedulerの設定・操作は変更していない。

### csh・ログ側の未完了項目

生成.csh 53,988本（約14.153 GB）、archive重複.list 53,990本（約14.184 GB）と、統合済み1CE成功stdout約26,941本（約7.06 GB）が追加整理候補。A76/1CEのSUMS tar SHA256はrunmetaと一致し、csh/listのbasenameとサイズの全収録を確認済み。1CEの失敗49件、代表成功ログ、A76全ログ、未統合2AFの記録は保持対象。Condor照会で元workspaceにactive workerなしを確認済み。

候補一覧の保存処理は自動承認review timeoutで初回・1回再試行とも未実行だったため、csh/list/logの削除は行っていない。既存[承認障害時手順](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/docs/ai/skills/approval-recovery.md)の「長い再試行ループを続けない」に従い停止し、ユーザーへ手動承認への一時切替を依頼済み。ROOT整理は独立して完了した。

[機械可読summary](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261109_001/summary.json)／[削除後検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261109_001/postcheck.json)／[削除一覧・検証根拠・実行journalの圧縮記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/storage_cleanup_20261109_001/audit_records.tar.gz)。削除計画SHA256: `fd4100b7be376e699f322abf997c79ed3f7b3fad567eb31a036ec7ee9c40df4c`。git commit/pushはしていない。


## 2026-11-11 続行

ユーザーが手動承認へ切り替え、未実施の.csh/list/成功ログ整理を再指示した。134,919 files /35.399 GBを対象にtmux起動済み。7.7 GeVの追加解析は保留。[11/11の記録](../../11/11/storage_cleanup_logs_20261111.md)。


## 2026-11-11 残件整理完了

手動承認で.csh・入力list・検証済み成功stdoutの134,919ファイルを削除し、35.399 GBを追加解放。quota87.21%、未確定使用量控除後の空き330.257 GB。保護ファイル・独立検査PASS。[完了記録](../../11/11/storage_cleanup_logs_20261111.md)。
