# anaFemtoLambda_4He_tree 本番投入

**最新（2026-11-18）**：全4,058 jobs終了、正常footer付き回収ROOT4,024、未回収34。最終2869も正常回収。ユーザーの追加依頼で[近接runを20群へhaddする処理](../18/tree_run_period_hadd.md)をtmuxで起動。初段終了と全件正常完了を区別し、元の投入時点の記録は以下に保持する。

ユーザーが2026-11-17に新tree解析の実行・job投入を明示指示。旧2AF削除後のquotaを確認し、compact初段tree＋108QAを作る。未較正81runは保留、SE/MEは保存ROOTを読む後段で行う。

状態：**本番4,058 jobsの全投入が完了**。SUMS ID **`9037A6F05FAD08FC8E9BD9286C8A38A5`**。2026-11-17T15:24:47.853130+09:00時点でCondor実行中1,413・待機2,469・Hold0、解析は実行中。73,054入力coverage、保存資料5種SHA、完了通知processの稼働を独立自動監査しPASS。

- 入力：2021 P24iy、較正済み388run・73,054 unique files。frozen list SHA256 `b469b4a926dd4781256d3e289c08655c53dd71f41f343411967df342070fd9c5`。
- mainconf：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml`。`storePrePidControl:false`、25files/subjob以下、maxEvents=-1。
- XML：`job/joblist/joblist_auau3p85_anaFemtoLambda_4He_tree_treeWrite20261117.xml`。
- 出力：`rootfile/auau3p85_anaFemtoLambda_4He_tree/production20261117/`。treeとQAは同一ROOT。元producerごと保持し、自動haddしない。
- generator概算2,923に対し、実SUMS balancing（min18/max25）は4,058jobs。完成後に全入力対応を自動検証する。
- 作業証跡：親の`rootfile/femto_lambda_tree_production_20261117_001/`。build/source比較・smoke closure・固定input・launch設定とログを保持。
- 初段QA108個・後段pair1186個・EP総和32比較は全bin/content/error/entries/flow一致。He4 SE13/ME403、EP1 ME459。
- 初段規模の参考外挿433.35 GiBは旧EP8共通成功4023groupのevent数を基準とする。今回全73,054入力の実Entriesではない。現user soft余裕約2.25 TBを確認し、1.2 TB以上を投入直前guardとした。
- 以前に承認済みの完了メールをminimal内容で設定する。通知は全job終了と出力/完了ログを確認して一度送る。投入・再投入・haddを通知scriptに行わせない。

- [x] batch build、fresh smoke/closure、入力/較正388一致、quota、XMLを確認。
- [x] 一回限定のtmuxから全4,058job投入を完了し、SUMS ID・実job数・queue稼働を確認。
- [x] runmeta/configlog/joblistlog/runtime/SUMS tarと入力分割coverageを検証し、完了通知の稼働を確認。4058contiguous .list、73054入力、maxFiles=19、保存資料SHA照合PASS。
- [ ] job完了後、全ROOTのtree/metadata/coverage・実容量を確認し、後段再解析へ進む。
- [ ] [Δφ*–Δηを再作成し、旧全統計PDFとの傾向一致を確認](old2AF_retirement_and_merging_recheck.md)。今回の投入前smokeだけでこの比較を完了扱いにしない。

commit/push・旧ROOT再削除・未較正run較正は今回行わない。


## 起動済みの継続処理と保持記録

- tmux `he4_tree_prod_20261117_001`：一回限りのsubmitを継続。入力coverage/runmeta照合後に完了メール監視を自動起動する。
- tmux `he4_tree_prod_audit_20261117_001`：全投入後の4,058contiguous .list、73,054入力coverage、runmeta/configlog/joblistlog/runtime/SUMS tarのSHA、completion-only監視の実PID/設定を確認し、`post_submission_audit.json`へ保存する。監査はqueue変更・hadd・メール送信を行わない。
- 実worker4057と固定ZIPの13critical項目（source/config/lib/dictionary/EP較正）のSHA一致を独立確認した。
- 最新確認時soft余裕は2.237 TB。quota不足や通知済みを推測しない。

[小さい検証・投入証跡](../../../../verification/tree_production_20261117_001/README.md)。原本live状態は親`rootfile/femto_lambda_tree_production_20261117_001/tree/launch_state.json`と`post_submission_audit.json`。実施中のSUMS IDがあるため、同じXMLを再投入しない。


## 全投入・保存記録・通知監視の確認完了

2026-11-17T15:24:47.853130+09:00に全投入後のphase=`monitor_started`、post-submission audit=`passed`を確認。通知PID 2819679 が実行中で、全job終了後にminimalメールを一度送る。現在のメール送信記録はFalse（投入完了と解析完了を区別）。確定runmetaは`job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_tree_treeWrite20261117_9037A6F05FAD08FC8E9BD9286C8A38A5.json`。全4,058入力listと実workerを保持し、native helperの再計算でもcoverage・分割順序が一致した。

現在ROOT回収163本、soft quota余裕は約2.141TB。新treeの全ROOT schema/entry/metadata/coverageと実サイズの最終評価はjob完了後のTodo。初段treeの自動haddは行わない。


### 実farm出力の再読確認

回収されたsubjob3405のROOT（30,749,291 bytes、保存118,240event）を元projectから再読し、10,000保存eventのEP8 replayを実行してstatus0・completed10000でPASS。scratch/runtime生成→GPFS回収→project側Reader/共通pair計算まで動作確認した。PicoDst再読なし、元ROOTのsize/mtime/inode不変。これは1本のsample検証で、全4,058出力の最終検証ではない。証跡は`farm_probe_result.json`。


## 2026-11-18進捗：未完了・未回収出力あり

[最新状況](../18/tree_job_status_20261118.md)。10:18 JST時点でRunning3、終了記録4055に対し回収ROOT4021。差分34は入力PicoDst読取り失敗等で未回収。全件正常完了とは扱わず、job取消・再投入はこの照会で行わない。
