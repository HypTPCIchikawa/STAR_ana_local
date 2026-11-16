# 3.85 GeVヒストのみでの再投入

状態：再投入開始済み（2026-11-04T15:41:57.536321+09:00確認）。SUMSはcatalog照会・入力分割を完了し、サブjob生成・投入を継続中。新IDは`1CE33994716F7896423A80AB994BAB14`（生成ログで確認。全投入終了後にrunmetaを保存）。2021 P24iyの149,817保存場所から重複を除いた80,982 PicoDst、予定26,994 subjobs。全件投入完了・全件解析完了とは区別する。

- [x] 容量不足を確認し、ユーザーが「候補保存なし・ヒストのみ」を明示選択。
- [x] candidateStore行だけを除いた専用mainconfを作成。jointSparse:falseを維持。
- [x] batch-matched SL24y buildと4,027イベント実行が正常終了（選択3,858）。
- [x] 1,786denseヒスト、59,146,854cellsの内容・誤差・Entries等が候補保存ありの参照と完全一致。
- [x] He4 signal SE=9/ME=251を確認。候補TTree/layout/modeおよび8D/metadata不在。ROOTは3,538,869bytes。
- [x] 2021 P24iy catalog XML、589ファイルのfingerprint、投入前quota budget、一回限定launcherを検証して開始。
- [x] 全件投入後に、今回のrunmetaだけを用いる最小メール通知へ移行するlauncherを起動。通知の16テストPASS。
- [x] 記録と日報を更新。小さな設定・検証結果・マクロは個人snapshotへ保存。

mainconf：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`。PID/KF/TOF/混合条件と通常ヒストを維持し、候補保存だけを停止。既存2AFの候補ROOTは保持し、新しいIDと混合しない。候補cacheがないので、後の再混合にはPicoDst再処理が必要。

投入直前soft quota空き372.674GB（blockInDoubt控除後）、要求250GB以上。新規ROOT約150–180GB、log/err/job約30GBを見込む。自動hadd・自動再投入なし。

永続実行：starsub03のtmux `he4_histonly_20261104_001`。主状態は`/tmp/yichikawa/hist_only_relaunch_20261104_001/launch_state.json`、提出ログは同directoryの`submit.log`。GPFSへの状態コピー失敗が主状態を失わせない構成。全件投入終了とrunmeta一致を確認後、`--email-content minimal --state-dir /tmp/...`で解析終了・要確認メールを既定宛先へ一度送る。現在は`submitting`であり、メール監視は全投入後に起動する。

[再現手順・運用記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/hist_only_relaunch_20261104_001/README.md)／[少統計の完全一致検証](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/hist_only_relaunch_20261104_001/closure/validation.json)。Git commit/pushなし。

稼働確認：新IDをUserLogで指定したCondorキューの状態は{"2": 211, "1": 189}（1=待機、2=実行中）。先行subjob26993の完成cshにも新mainconfが指定されていることを確認。独立監査でも設定・一回限定投入・通知への移行に問題なし。


## 2026-11-05 16:37 JSTの更新

全件投入は完了。26,987件に終了記録（exit_0=26,979、removed=8）、7件が実行中。元launcherのreport ID読取エラーによりメール監視が起動していなかったため、保存済みrunmetaを検証し、今回IDを指定して最小通知監視を直接起動・稼働確認した。元の失敗記録を保持し、復旧記録を別途保存。[最新状態と復旧記録](../05/job_status_notification_recovery_20261105.md)。7.7 GeV準備は現行job終了後まで保留。


## 2026-11-06の更新

全件終了・queue0と、11月5日18:56 JSTの終了メール受理を確認。回収26,945 ROOT（49件除外）を全件metadata/key検査し、実測event数が旧04BF比約−0.045%の小差であることを確認。ユーザー指示でtmux `he4_hadd_1ce_20261106`に統合処理を起動し、hadd終了後のメール通知を設定した。hadd完了とはしていない。[最新Todo](../06/hadd_1ce_20261106.md)。
