# 保存treeの期間内連続ME・pool100本番

状態：**全20job正常完了・ROOT20本回収済み。2026-11-20 21:28〜21:34 JSTに完了ログ・出力件数/設定を照合してPASS。** 最終job終了は11/19 23:12:15 JST、完了メールは23:20:55 JSTにlocal sendmail受付済み。SUMS ID `2D07F680577123CABC16FDD9159CFD9E`。

投入時の記録：2026-11-19 16:10 JSTに20jobと保存記録・coverage・queue・完了通知監視を確認。直後はRunning11・Idle9・Held0。ユーザーが期間内でpoolを継続し、mix対象event上限20→100への変更・動作確認・本番投入を明示依頼。

- 入力：検証済み20期間ROOT、計4,024 producer・保存1,819,447,887 events。1期間1job、期間間のmixingは行わない。
- mainconf：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_periodMix100.yaml`。
- `resetPoolAtSegmentBoundary:false`、`requireOriginalMixing:false`。元producer境界でpoolを初期化せず、各job開始時に空のpoolを用意する。
- 専用mixing YAMLで`bufferSize:100`・`maxMixEvents:100`。cent9×vertex10×EP8区分、bufferAll・双方向を維持。各poolで最大100 admitted eventsを使用し、常に100個の相手候補が存在することは保証しない。
- 保存済みcompact候補だけを再読し、PicoDst reader・Λ再構成は起動しない。
- 作業証跡：`rootfile/femto_lambda_period_mix100_20261119_001/`。

- [x] 既存codeのpool制御と二つのdepth上限を確認し、専用4 YAMLを用意。
- [x] batch整合build・設定読取を確認。初回replayで見つかった`vector<float>`初回bindエラーを、project LinkDefへのfloat/double/int/unsigned int辞書明示登録で修正。再build後、補助ACLiCを挟まない通常runnerのfresh processで1event読取りに成功。
- [x] 同一保存treeの先頭400,000eventsをB=継続100、A=reset100、C=継続20で比較。3segment、境界1/182434/357248を実通過し、SE607hist・flow込み28,242,428bin比較でcontent/error/entries/axes完全一致。MEはB602,979/A571,563/C123,864 pairs、cap除外0、両方向非空。
- [x] Λ–⁴He signal MEはB116,183/A109,474/C24,000。継続20→100で約4.84倍、reset100→継続100で約6.13%増。400k試験の値で、全統計の増加率を保証しない。
- [x] quota、入力20件coverage、worker staging/メモリ/時間、保存runtimeのSHAを確認。最新preflight空き1,229,916,471,296 bytes、最低500GBをguard。native `request_memory=4096` MB、`request_disk=65536000` KB。各workerは期間ROOTをscratchへstageして読む。
- [x] mainconf指定のjoblistを直前生成し、tmux `he4_period_mix100_20261119_001`から1回だけ20job投入。SUMS ID・runmeta/configlog/実runtime ZIP・20 native入力list・live queueを確認、minimal完了通知PID1776689を起動。
- [x] 新depth100のEP区分総和も既存PlaneSums helperで確認。He4/He3×4mass gate×SE/MEの16 pair（32hist）、flow込みのcontent/error/Entries/axis差0。
- [x] 全20jobでrequested=completed=available=各periodの予定保存event数を確認。合計1,819,447,887events、Maker completed・正常native exit、StChainのWarn/EOF/Err/Fatal counterはすべて0、打切りなし。
- [x] 回収ROOT20本、計11,491,533,504 bytesを確認。全ROOTがcompleted/replay、false/false・両depth100・bufferAll/bidirectional。eventLedger件数と先頭/末尾のevent・元producer範囲が一致し、samplerのcap除外0。確認用に全eventを再解析したわけではない。
- [ ] 20期間出力を用い、SE/ME・CF・mergingを比較する。高k*でCFが1へ近づくかの物理判断は今後の作業。

[未回収34件は保留](unrecovered34_cause_investigation.md)。未較正81runも保留を維持。原本4024ROOTと20期間ROOTは保持し、既存設定を上書きしない。commit/push・自動再投入は行わない。

## 実測・本番の参照先

- B400kのevent loop：wall108.744秒/CPU98.48秒、wrapper全体wall2分08.57秒、peak RSS544,092KiB。C20のloop CPU36.82秒に対し約2.67倍。最初の1periodのprefixによる実測であり、全runの上限・ETAではない。
- 本番joblist：`job/joblist/joblist_auau3p85_anaFemtoLambda_4He_tree_periodMix100_periodMix10020261119.xml`。XMLの`minFilesPerProcess=maxFilesPerProcess=1`と`filesPerHour=0.025`を明示。後者はscheduler向け40時間/periodの見積で、処理時間の実測ではない。
- runmeta：`job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_tree_periodMix100_periodMix10020261119_2D07F680577123CABC16FDD9159CFD9E.json`。同IDのconfiglog/joblistlog/runtime_bundle/sums_artifactsを保持。
- 出力：`rootfile/auau3p85_anaFemtoLambda_4He_tree_periodMix100/production20261119/`。
- live状態：studyの`launch_state.json`・`input_coverage.json`・`native_resources.json`・`notify/{settings,status,sent}.json`。通知は全20件のterminalと出力/解析完了markerを確認して一度だけ送る。自動hadd・再投入は行わない。
- 16:16 JSTの独立queue照会では20件すべてRunning・Held0、live RequestMemory/RequestDiskも指定値に一致。投入完了と解析完了は区別する。
- 初回dictionary失敗の3pilotは`initial_dictionary_failure/`へ指紋付きで退避。修正後の`pilot_status.json`・`pilot_comparison.json`とは区別する。原本ROOTは変更していない。

[個人submodule内の検証・実投入証跡](../../../../verification/period_mix100_20261119_001/README.md)。ROOT・binary・bulk logの共有mainへのcommit/pushはしていない。

## 本番完了確認

- 完了通知記録：`notify/status.json`は20 submitted/20 terminal/queue0、exit0が20、正常header20、missing/badHeader/analysisCompletionUnconfirmedすべて0。`launch_state.json`は`notification_completed`。
- 各jobのevent-loop wall時間4.02〜6.96h、CPU3.99〜6.87h、計97.24 CPUh。最後はgroup00/subjob0が11/19 23:12:15 JSTに終了。
- Λ–⁴He signalのraw pair Entries合計：SE2,221,923、ME459,780,318。Λ–³He signal：SE10,842,199、ME2,141,220,417。規格化・背景差引き後のCFではない。
- ROOT監査scope：全20fileのmetadata/counters、設定、signal hist・mixing QAとledger先頭/末尾40rows。実解析が全保存eventsを処理したことは、独立した全20完了ログと予定数の照合でも確認した。
- 完了監査：study内の`completion_verified_20261120.json`、`completion_log_audit_20261120.json`、`output_root_check_20261120.json`。固定出力一覧は`output_roots_20261120.list`。
- [個人repo内の完了証跡](../../../../verification/period_mix100_20261119_001/completion_20261120/)にもSHA照合して保存した。今回の確認で出力hadd・新たな解析job・メール再送は行っていない。
