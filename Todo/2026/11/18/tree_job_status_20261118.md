# 新tree本番jobの状況（2026-11-18）

**最終：全4,058 jobs終了、queue0、回収ROOT4,024本。未回収34件は前回と同じ入力エラー等。** 最後の2869は正常回収。17:39 JSTに要確認メール送信記録あり。過去時点の記録は下記へ保持。[run時期別haddと欠損監査](tree_run_period_hadd.md)へ継続。

10:18 JST確認。SUMS `9037A6F05FAD08FC8E9BD9286C8A38A5`、予定4058jobs。**未完了**。

- Condor終了記録4055件（外側exit0）、Running3件、待機/Holdなし。
- Runningはsubjob1476/2869/3788。18〜19時間経過しているが、子ログはscratch内にbufferされるため共有stdoutだけでは進捗・停止を判定できず、ETA未確定。
- 回収ROOT4021本、すべて非空。logical 471,108,138,019 bytes。全内容/schemaの検査済みという意味ではない。
- 終了記録と回収との差34件には解析完了markerがなく、32件がunreadable PicoDst、1件がempty PicoDst、1件がchain Make=2途中終了（requested502094/completed501031/status4）。
- SUMS wrapperのexit0だけで解析成功としない。未回収34件は完了扱いにせず、入力データ/転送/途中読取の確認と再処理要否の判断が必要。
- 完了通知processは稼働中。全job終了前のためsent記録なし。最終通知は不足出力があれば要確認になる。

今回の照会では読取りと記録だけを行い、job取消・再投入・コード変更・メール送信は行っていない。

[投入Todo](../17/tree_production_20261117.md) / [今回の明細](../../../../verification/tree_production_20261117_001/status_20261118_1018.json)


## 14:36 JSTの再確認

前回残り3件のうち1476は10:55、3788は12:36にrequested=completed/status0で正常終了し、ROOTを回収。残る2869（Condor794176.88、spool1801）はRunning。共有stdout376B/stderr0Bのままで、内部の進捗・停止は判定できない。終了時刻は未確定。

終了記録4057に対し回収ROOT4023本、logical471,349,869,260 bytes。前回判明した未回収34件は変わらず、正常全件完了ではない。完了通知processは稼働しsent記録なし。取消・再投入・メール送信はこの照会で行っていない。

[14:36の明細](../../../../verification/tree_production_20261117_001/status_20261118_1436.json)。


## 16:59 JSTの再照会と最終jobの内部確認制限

終了記録4057/4058、queue1（{'2': 1}）、回収ROOT4023本。通知のsent記録はFalse。最後に確認したRunningは2869/Condor794176.88/spool1801で、約26時間経過。共有stdout376B/stderr0Bのまま。

直接SSHは既知host keyがなくstrict checkingで停止し、検証を無効化していない。正式Condor経路の内部ログ読取も`SCHEDD_ENABLE_SSH_TO_JOB`で拒否された。サーバー側policyを迂回せず、解析eventの進捗・停止は未判定、ETA未確定とする。CPU/leaseの更新だけで解析進行を断定しない。jobの取消・再投入・解析コード変更・手動メール送信なし。

[この時点の記録](../../../../verification/tree_production_20261117_001/status_20261118_1659.json)。

## 最終確認：メールの「正常完了を確認できないjob」

通知の17:32 JST snapshotでは全4058 terminal、queue0。2869はrequested=completed550829、reconstructed543693、status0でROOT回収済み。未回収34は32 unreadable PicoDst・1 empty PicoDst・1 Make途中終了。4058件の外側exit0と4024件の解析正常完了を区別した警告である。minimal通知は17:39 JST local sendmail受付済み。

今回34のworker入力groupは前回EP1/EP8と全て一致し、33件は同じ問題Picoまで特定できた。Make途中終了1件は原因Pico不明のまま。612/73054入力Pico（0.837736%）が未回収workerに属し、うち310本の実測Entriesは合計8,594,396、302本は不明。欠損入力割合の下限0.464210%であり、accepted/pair割合ではない。run別最大36.69%の偏りがあるため、統計的に無視できるとは断定しない。再投入は未実施。詳細は[未回収入力監査](../../../../analysisnote/20261118/tree_missing_inputs_20261118.md)。
