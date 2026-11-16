# 3.85 GeV job状態と通知監視の復旧

確認：2026-11-05 16:37 JST。SUMS ID `1CE33994716F7896423A80AB994BAB14`。

全26,994件の投入とlistの連番0–26993を確認した。終了記録26,987件（exit_0=26,979、removed=8）、実行中7件。残りはsubjob 91、722、813、926、1537、1953、6390で、確認時の実行経過は約21.6–22.1時間。進捗が正常か停止しているかは未判定であり、キャンセル・再投入は行っていない。

ROOT名は26,938本、未回収56本。そのうち7本は実行中分、49本は終了記録があるが出力がない分。標本3本の先頭はROOT magicだが、全内容/event数の検査はまだ行っていない。exit_0だけで全件正常と判断しない。

## メール監視の不備と復旧

元launcherは2026-11-04 19:28 JSTにsubmit.shのreturn 0を受け取った後、SUMSのreportファイル名からIDを読み取れず停止していた。runmetaではreportPath/sessionPathがnull。一方、submit.sh側はprocess IDから今回IDを取得して、config/joblist snapshotとrunmetaを保存済みだった。元の失敗記録は`/tmp/yichikawa/hist_only_relaunch_20261104_001/launch_state.json`に保持する。

runmetaのID/mainconf、全list coverage、Condorイベントとの一致を確認し、既存の通知スクリプトを今回のrunmeta指定で直接起動した。通知スクリプト自体は検証時のfingerprintと同一で、解析のソース・設定・libraryは変更していない。二重監視/送信記録の有無を確認してから起動した。

復旧後のtmuxは`he4_histonly_notify_20261105`、PID304822。16:37:18 JSTのstatusにsubmitted=26994、terminal=26987、queued=7を確認。宛先は`yudai.ichikawa.d3@tohoku.ac.jp`、`--email-content minimal`、5分間隔/settle5分。監視が現在稼働中であり、全job終了後にROOTの回収状態も確認し、完了または要確認メールを一度送る。復旧確認時点で送信はまだ行われていない。8件のremovedと欠損があるため、全件正常完了とは扱わない。

主監視directoryは`/tmp/yichikawa/hist_only_relaunch_20261104_001/notify_1CE33994716F7896423A80AB994BAB14`。復旧記録は[notification_recovery_20261105.json](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/hist_only_relaunch_20261104_001/notification_recovery_20261105.json)。過去の「全投入後に監視へ自動移行」という設定は、元launcherの失敗により実行されなかった。今回の直接起動により運用を復旧した。

- [x] 実行中jobと終了・出力coverageを確認。
- [x] 正しいrunmetaから最小メール監視を起動し、稼働statusを確認。
- [x] 全job終了・終了内訳・ROOT回収を確認（11月6日）。
- [ ] 回収ROOTの全内容をhaddのcheckpointで検証する。
- [ ] 今後launcherを再利用する前に、report出力がない場合のID取得をrunmeta等との一致検証付きで修正する。

自動hadd・追加job投入・キャンセル・ファイル削除・Git commit/pushは行っていない。7.7 GeVの準備はユーザー指示どおり現行job終了後まで保留を維持する。


## 11月6日の更新

全26,994件が終了しqueue0、回収26,945/未回収49を確認。復旧した監視から11月5日18:56 JSTに終了メールが送信サービスへ受理された。ユーザー指示で部分haddをtmux起動し、hadd終了後のメールも設定。[haddの最新状態](../06/hadd_1ce_20261106.md)。
