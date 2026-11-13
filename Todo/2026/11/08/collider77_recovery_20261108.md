# 7.7 GeV未処理分の再投入

状態：11月8日15:56:31 JST、tmux `auau7p7_recovery_20261108`で分割再投入を開始。最初の1200入力／200 jobsは投入完了・監視中。残りを全件投入済みとは扱わない。進行中の2090 ROOT部分hadd・QAには干渉しない。

## 対象

- 未開始でremovedとなった2954元jobs／17725 URI。
- 入力不良10元jobs／60 URIは再取得・検査した有効入力を再処理し、未解決入力を別途記録する。
- bad-run-only31元jobs／186 URIは再投入しない。既存のbad-run cutを保持する。
- 既存ROOTに対応する12540 URIとの重複を禁止する。

回復候補17785 URI＋既存12540 URI＋bad-run-only186 URI＝元30511 URIを既に照合済み。job再分割により新subjob数は元2964件と必ずしも一致しない。

## 容量・投入方式

11月8日15:34 JSTの空きは115.95 GB（soft quotaからusageとinDoubtを差引）。再投入40 GB（ROOT30、metadata/log7、その他3）、進行中hadd追加40 GB、安全余裕20 GBを予算とする。次batchのrunning分についても1 job25 MBを予約し、約200 jobsで追加約5 GBを条件に含める。これはモデル上の保守的予算であり、絶対上限ではない。

原則最大1200入力／約200 jobsずつ、前batchの全件terminal・queue0・producer/ROOT検査後に次を投入する。未投入分はローカルmanifestで待機し、Condor hold/release/removeは使わない。quota不足時も新規投入を待つだけで、実行中・旧jobを操作しない。自然発生のheld状態は期限前に検知して要確認とし、強制releaseしない。

未開始17725 URIは14×1200＋925の15batch。入力不良60 URIの検査完了後、有効分を追加batchにする。全batchと入力分類が確定するまで完了通知を送らず、最後に一回だけ最小メールを送る。自動の失敗job再試行ループを作らない。

## 作業

- [x] 現quotaと進行中haddの状態を確認。元部分統合は15:31時点で300入力の検証完了、継続中。
- [x] 独立した容量modelで再投入ROOT約15–15.5 GB＋metadataを見積もり、40 GB予約を採用。
- [x] 既存SL22b build・source/config fingerprintを照合。物理選択・PID・TOF・8D/cache無効を変更していない。
- [x] 未開始15batchのXML・unique ZIP packageを作成し、実SUMS simulationで16入力→8＋8、request_disk20 GB、URI coverage一致を検査した。simulationはdispatchなし。
- [x] 入力不良60 URIをbounded cacheで再取得し、ROOT5でheader/schema/全basketを検査した。51本が再処理可能、3本は正常な0-event入力、6本は未解決。健康な別replica1本を回復した。
- [x] 分割controllerの15 mock casesとmetadata adapterの6判定cases・comment除去・bad-run-only統合testを検証し、tmuxで一回限定起動した。testは実job/メール操作なし。
- [x] 最初のSUMS ID`482C1A5B79B005C132DD7BEEE5AFADEC`、200件の全投入・URI照合を確認した。16:00:47 JSTにrunning124、idle74、terminal2（exit_0）、held0。
- [ ] 全batch終了後に実処理coverage、正常ROOT・producer完了、未解決入力を記録してメール通知する。
- [x] Todo/日報と運用・検証記録121ファイル計3,735,960 bytesを個人submoduleへ保存した。ROOT・大ログ・libraryはコピーせず、commit/push・元ROOT削除はしていない。

## 終了判定

SUMS外側wrapperは内側解析失敗をexit_0で隠す場合があるため、Condor終了だけを正常判定にしない。既存のproducer footer requested/completed/statusとROOT headerをbatchごとに確認し、未完了を記録する。外側SUMS shellを生成後に書き換える方式は採用せず、元template・解析codeを保持する。

ROOT5の直接root://読み取りはplugin環境のため初回検査が失敗した。この結果を60入力のデータ破損とは扱わず、1本ずつ自分専用の一時cacheへstageして検査する方式へ切り替えた。cache予算2 GB、元データは削除しない。

作業領域：独立workspaceの`rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/recovery_20261108_001/`。

## 起動記録と入力検証

起動前15:55:56 JSTに空き115.54 GB、次200 jobs用の必要予約105 GB（再投入40＋hadd40＋安全20＋writer5）を確認。source/config/15 library・初期15 XMLと入力list・検査マクロをhash固定し、haddの進捗状態や未確定入力検査のrolling結果は固定hashへ含めない。

主状態は`/tmp/yichikawa/recovery77_20261108_001/status.json`、GPFS mirrorは新作業領域の`controller/status.json`。controller PID1429849。初期15 batchesを順番に処理し、最後に入力不良60本の確定planを読み、全60本が有効追加batchまたは明示除外に分類されるまで最終通知を出さない。各batchのheader/footerとROOT metadataを照合し、外側exit_0だけでは正常としない。

51本の実測PicoDst Entriesは708,816。うち`st_physics_22110045_raw_2000002.picoDst.root`はcanonical経路の不良copyを避け、qgp008:36607上の同一LFNの健康なcopyを検証した。hostname形式で実xrdcp＋ROOT5全7172 basketsの読取りを通過、17,663 entries／600,655,270 bytesを確認。残る6本は未解決（Zombie等5、basket破損1）、3本は20 branchesを持つ正常な0-entry ROOT。catalog件数を実測件数として使わない。

追加batchではlogical input URI/LFNの同一性を保ち、この1本だけworkerのUSERCODEで検証済みPFNから再stageする。元SUMS wrapperによるcanonical staging後、readability guard前に、exact LFNだけをprivate FILELISTへ置換する。元template・初期15 XML・実行中controllerは変更していない。追加batchはpending planへ登録済みだが、まだ投入済みとは扱わない。

## 最終準備完了・稼働確認

- [x] 追加51 canonical URI＋除外9 URI（正常empty3／未解決6）で元60 URIを完全に分類。全回復候補17785 URIの重複・欠落なしを、稼働中とは別の一時controller状態で検証した。
- [x] 追加v2 XMLはplain xrdcp -fを最大3回のmanual loopで実行し、各回にstatus・600655270 bytes・ROOT magicを検査。workerのstarver後に未知の--retry optionへ依存しない。旧v1のcode/XML/planは履歴として保存。
- [x] 実C-shellの8 mock casesがPASS。転送2回失敗→3回目成功、連続失敗、size/magic不良、非対象、重複、元入力不変を確認。実ROOTデータ転送はtestで行わない。
- [x] 実SUMS simulationで追加51入力→8 jobs（7,7,7,6,6,6,6,6）、全request_disk20000000 KB、全URIが一回、dispatch/queue/eventlogなし。追加batch上限9 jobs。
- [x] 初回batch監視中・pending未読を確認し、16:17:37 JSTにprimary planをv2へatomic切替。稼働中controller/source/config、初期15 batchesは変更していない。

primary plan SHA256：`47e3219232ee6e7cf6fdf05cd3d2c6a378620bc792c152fa6859fa4fe5f32bf1`。controller SHA256：`7b967c9baf4c76f78dcd48b248c49281fc05eb608a2c96c317990112736748b3`。

16:20 JSTの主状態はmonitoring。初回200件のうち37 terminal（exit_0）、163 running、held0を確認した。回収済みsubjob30 ROOTは27,952 input／3,068 reconstructed／Lambda140、hist1786・ledger schema・producer footer一致をROOT5で検証済み。元のROOTは不変。FullField18入力の3 groupsはrun22062034–22062036が全て既存BadRun.hに入り、全入力処理済み・選択0・status4となる既知のゼロ選択と確認し、故障として再試行しない。

同時刻、元2090 ROOTのhaddは1300入力／13 checkpointsまで完了し、処理継続中。今回再投入の新IDはその入力listへ混ぜない。全再投入batch終了とproducer/metadata検査後に一回だけメール通知し、未解決入力や出力があれば要確認とする。

[再投入の説明](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/recovery_20261108_001/README.md)／[追加v2 plan](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/recovery_20261108_001/pending_validated_inputs_plan.json)／[主状態mirror](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/recovery_20261108_001/controller/status.json)。

16:24:53 JSTの初回batchは200 submitted、42 terminal（exit_0）、158 running、held0。全batchや追加51入力が既に終わったとは扱わない。[個人submoduleの保存snapshot](../../../../workspaces/auau7p7_2021_P22ib/recovery_20261108_001/snapshot_manifest.json)は保存時点の履歴で、live状態は上記主状態から確認する。

[終了通知の調査・部分hadd](collider77_hadd_qa_20261108.md)／[元manifestと分類](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/hadd_qa_20261108_001/completion/recovery_manifests.md)。


## 状況確認 — 2026-11-09T14:44:15.404000+09:00

旧ID `36B7E875A19554DB711364BE0B072DA5` の部分hadd・出力検証・QA作成の完了を確認。回収2090/5085 jobsを統合し、2026-11-08 16:53 JSTにQAと通知を完了した。22 checkpoints、1786 hist、全payload readbackと段階加算、event/TTree件数の検証がPASS。入力168,680,371、選択後44,203,571、Λ候補13,729,465。再投入分は未統合。

[QA PDF（10ページ）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/auau7p7_qa_partial_36B7E875A19554DB711364BE0B072DA5.pdf)。RFF・run22090006–22121018に偏った部分統計。選択後Vzの範囲外/flow0、mixing Vz clamp0。4He候補4227（低R4108、高R＋TOF m²119）。PID再較正・純度・CF評価は行っていない。[Λ質量・アクセプタンス](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/preview_20261109/qa_page_03.png)／[4He β–η・TOF条件比較](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/preview_20261109/qa_page_07.png)。βは候補の|p|/Eであり、TOF実測βとは区別する。

14:33:56 JSTの独立Condor確認では再投入controller稼働、初期batch1–9の1800 jobs終了、metadata監査でROOT1575本、既知bad-run全入力/選択0が197 jobs、要確認28 jobs。28件はproducer完了footerなし、詳細原因は未確定。batch10は200 submitted、67終了、133 running、idle/held0。未投入5776入力（初期5725＋追加51）。終了メールは全batch後で未送信。再投入ROOTの全payload検証は後続hadd時に行う。

今回のquota実測：使用2.486998 TB / soft 2.588045 TB（96.10%）、hard 3.137801 TB、inDoubt 6.595 GB、inDoubt差引き空き94.452 GB。部分hadd終了でhadd予約0、再投入40＋安全20＋次200jobs5＝65 GBの投入条件を満たす。各batch前に再照会し、容量不足時はscheduler投入前にローカル待機する。holdは使用しない。

[確認記録JSON](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/auau7p7_2021_P22ib/status_20261109_001/status_check.json)。解析source/config/library・controller、schedulerは変更していない。ROOT/ログ削除、git commit/pushは行っていない。追加stderr診断の一部は自動承認review timeoutで未実行であり、28件の原因を確定したとは扱わない。


追記：保存した最新controller snapshotでは 2026-11-09T14:43:15.591813+09:00 に第10batchは終了82、running118、idle/held0。上記14:33の数字から進行している。quotaは 2026-11-09T14:44:15.404000+09:00 照会で、inDoubt 6.595 GBを差し引き、空き94.452 GB。必要65 GBを満たす。


## 2026-11-11 ユーザー指示による保留

7.7 GeVの追加調査・解析・再投入・新しいhadd/QAを一旦保留する。既存controllerは11/10 11:57 JSTにcompleted_with_issuesで終了済みと読取り確認し、今回scheduler操作は行っていない。7.7 GeVの全成果物と記録を保持し、元3.85 GeV作業領域の.csh/list/成功ログ整理だけを完了した。[今回の作業記録](../../11/11/storage_cleanup_logs_20261111.md)。
