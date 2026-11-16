# 7.7 GeV終了通知の調査・部分hadd・QA

状態：2090本のmetadata検査を通過し、11月8日15:17 JSTにtmuxで部分haddを開始した。内容検証・元ROOT整理後に10ページのQAを自動作成し、終了メールを送る。まだhadd・本番QAの完了ではない。

初回停止の履歴：追加検査が自動承認レビューで初回・再試行とも実行前timeout、通常sandboxもPython起動前に`bwrap: Creating new namespace failed: No space left on device`で失敗した。ユーザーが手動承認へ切り替えて再開し、その後自動モードに戻して作業を継続した。namespace障害を解析ROOTの容量不足と混同しない。

解析：`auau7p7_anaFemtoLambda_4He_KFParticle`、2021 P22ib。
SUMS：`36B7E875A19554DB711364BE0B072DA5`。
作業場所：`star-analyzer_auau7p7_2021_P22ib`。

## 確認済み

- [x] 通知主状態はnotification_completed、11月8日09:51 JSTにlocal mail serviceで受理。
- [x] Condor submitted/terminalは5085/5085、queue0。
- [x] exit_0=2131、removed=2954（58.09%）。全統計の正常終了ではない。
- [x] 通知時のROOT header確認済み2090本。ただし現時点の全ROOT内容・producer完了・実event coverageを検証したものではない。
- [x] 通知監視の設定・状態と、capacity guardの終了状態を読み取った。guard最終記録では空き126.47 GB、旧hadd予約0だったが、これは現在quotaの再測定ではない。

## 再開後の作業

- [x] Condor009削除block/RemoveReason、未回収ROOT、producer未完了の内訳を調査した。大量の番号を本文へ並べず機械可読listへ保存した。
- [x] 30511入力URIに対する生成listの重複・欠落を確認し、回収/除外した入力・run/field範囲を記録した。catalog eventsを実処理数として使わない。
- [x] producer completed、event counters、通常1786histのkey/class/schema、通常eventLedger/selectedLambdaのentry/branch schema、候補cache/8D不在を全回収ROOTで検査した。
- [x] 実quotaと旧3.85 GeV hadd完了・staging使用量を再確認し、7.7 GeV統合用の各階層・最終検証・一時領域を予算化した。
- [x] 既存checkpoint driverを再利用し、今回IDの検証済みROOTだけを対象とするtmux処理を起動した。部分統計・除外listを明記し、既存ROOTを上書きしない。
- [ ] 統合ROOTの全payload読み取り、hist加算、event/TTree数・schema、入力不変を確認する。検証済み元ROOTだけをsubjobs/<jobid>/へ整理し、削除しない。
- [ ] QA PDFを作成する：選択前後Vzとshape比較、centrality/refmult/TOFmatch、Λ質量とacceptance、核種PID dE/dx/m²、4Heのlow/highR各TOF条件別beta–eta/beta–Y/Pt–Y、SE/ME k*分布。部分統計とrun範囲を全体に明記する。
- [ ] QAでは巨大eventLedgerの全行をGetEntries(cut)等で走査せず、検証済みmanifest/counters、tree metadataとhistを用いる。新しいCF/purityの結論は出さない。
- [x] 今日のREADME/日報と11月6日の投入Todoを結果へ更新し、運用コード・調査記録33ファイル計4,927,860 bytesを個人submoduleへ保存した。本番QA図はまだ未生成。git commit/pushはしていない。

removedが多いため、全統計が必要なら原因修正と未処理入力の再投入が別途必要になる可能性がある。今回は再投入・job設定変更・データ削除を行っていない。hadd/QAも開始済みとは報告しない。

## 根拠・再開点

通知記録：`/tmp/yichikawa/auau7p7_production_20261106_001/notify_36B7E875A19554DB711364BE0B072DA5/status.json`。
Condor event log：`/tmp/yichikawa/sched36B7E875A19554DB711364BE0B072DA5.condor.log`。
runmeta：新workspaceの`job/run/runmeta/runmeta_auau7p7_anaFemtoLambda_4He_KFParticle_36B7E875A19554DB711364BE0B072DA5.json`。

既存手順の追加読み取り自体が未実行のため、その内容を読んだとは扱わない。[承認障害時の手順](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/docs/ai/skills/approval-recovery.md)に従い、同じ要求を再試行し続けず、手動承認による再開へ切り替える。

## 手動承認後の確認

14:55 JSTの最新quota空きは116,937,719,808 bytes（116.94 GB、inDoubt差引済み）。回収ROOTは2090本、論理10,898,719,299 bytes／割当11,161,042,944 bytes、最小3,176,950／最大6,144,941 bytes。旧3.85 GeV haddは11月7日00:40:55 JSTにmerged、その後00:45:44 JSTに元ROOT整理まで終了しており、旧haddの追加予約は不要。

removed2954件のCondor履歴は全件000→012→009で、実行開始001がない。hold理由は今回の`STAR_CAPACITY_GUARD_36B7E875A19554DB711364BE0B072DA5`で、容量保護のために待機jobをholdした後、SUMSの「held状態が2日を超える」PeriodicRemove条件に達したことを確認した。解析コードが実行中に落ちた2954件ではない。容量保護とhold期限の組合せを見落とした運用上の問題であり、正常終了や単なる統計不足として扱わない。再投入はまだ行っていない。

既存再統合・日報skillを読み取り、checkpoint driverとROOT5検証を今回ID専用設定で再利用する方針を確認した。元ROOTは削除しない。10月のfloat差例外や3.85 GeVのevent比較判定を7.7 GeVへ流用しない。

## 調査・検査の完了とhadd起動

- [x] Condor履歴・除外2995件の理由・全30511 URIの重複なしcoverageを確定。
- [x] 全2090 ROOTのmetadata検査PASS。input168,680,371、reconstructed44,203,571、Lambda13,729,465、hist1786、eventLedger/selectedLambda schema一致、候補cache/Sparseなし。
- [x] 新規統合40 GB＋安全余裕20 GBを予約。15:16:45 JSTの実空き116.75 GB（inDoubt差引）で起動前check PASS。
- [x] 15:17:23 JST、tmux `auau7p7_hadd_qa_20261108`で部分hadd・QA・通知wrapperを一回限定起動。
- [ ] 21＋1の全22 checkpoints、最終payload読み取り・hist加算・event/TTree数と元入力不変を検証し、ROOTを公開する。
- [x] QAマクロ準備・ROOT5 pilot検証PASS（10 pages、41 QA hist）。通知wrapperの成功／merge失敗／QA失敗／重複起動保護をmockで検証PASS、実mailは送らない。
- [ ] 本番統合ROOTからQA PDF/JSON/PNGを作成し、部分統計の最小メールを送る。まだ完了・送信済みではない。

回収12540入力は全てReversedFullField、run22090006–22121018（1012 runs）。未開始17725入力は早いrunに偏り、FullField18本も全て未処理。従って、今回のQAは後半run/RFFに偏った部分統計であり、全期間やFullFieldのQAにはならない。catalog eventsは実処理event数と区別する。

ROOTがないexit_0の41件は、31件がbad-run-onlyで全event選択拒否・reconstructed0、5件がstaged入力header不良、3件がempty入力、2件がbasket/zlib不良による早期EOF。外側batch wrapperのexit0はproducer完了を保証しない。bad-run-only31件を物理解析失敗として再投入しない。破損が元データかstaging由来かは未確定。

容量guardではrelease条件を満たした2837 snapshotsでも自分のholdをguardHeld0と誤認して再開しなかった。HoldReasonの完全一致判定が認識を妨げた。native user注記付き文字列はregressionで再現できるが、削除後の歴史的ClassAd値そのものは未取得であり、suffixを直接測定済みとはしない。修正版copyと16 mock testsを新しいrecovery_preparation内に作成しPASS、旧guardは保持し、適用・scheduler操作・再投入はしていない。長期hold期限対策として将来は容量に応じた限定投入が必要。

主状態：`/tmp/yichikawa/hadd_qa_77_20261108_001/operation_status.json`、GPFS側は[operation_status.json](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/hadd_qa_20261108_001/operation_status.json)。最終haddとQA両方が成功した後にだけ通知する。途中失敗は要確認通知とし、自動retry・自動再投入を行わない。

回復候補の入力manifestも作成した。未開始2954 jobs／17725 URIと、入力不良10 jobs／60 URIを合わせて2964 jobs／17785 URI。bad-run-only31 jobs／186 URIは別listとし、回収12540 URIとの和が元30511 URIに重複・欠落なく一致することを検証した。再投入は行っていない。[回復候補manifestの説明](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/hadd_qa_20261108_001/completion/recovery_manifests.md)。

続報：その後ユーザーが再投入を明示指示し、15:56 JSTに別tmuxで分割再投入を開始した。元の部分hadd入力は変更せず、新IDを混ぜない。入力不良60本は再検証により51本を追加対象、0-event3本と未解決6本を明示除外した。[再投入Todo](collider77_recovery_20261108.md)に以後の状態を記録する。

15:22:09 JSTに最初の100 ROOTの統合・全payload読み取り・加算検証が通り、1 checkpoint完了、part1へ進んだ。実hadd子プロセスの稼働も確認済み。wrapper PID1323905、driver1323987。全体22 checkpointsと本番QAは引き続き処理中。

[個人submodule内の保存記録](../../../../workspaces/auau7p7_2021_P22ib/hadd_qa_20261108_001/snapshot_manifest.json)。ROOT・巨大ログ・ライブラリはコピーしていない。

[調査・実行の説明](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/hadd_qa_20261108_001/README.md)／[metadata事前検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/hadd_qa_20261108_001/merge/preflight/full_summary.json)。


## 状況確認 — 2026-11-09T14:44:15.404000+09:00

旧ID `36B7E875A19554DB711364BE0B072DA5` の部分hadd・出力検証・QA作成の完了を確認。回収2090/5085 jobsを統合し、2026-11-08 16:53 JSTにQAと通知を完了した。22 checkpoints、1786 hist、全payload readbackと段階加算、event/TTree件数の検証がPASS。入力168,680,371、選択後44,203,571、Λ候補13,729,465。再投入分は未統合。

[QA PDF（10ページ）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/auau7p7_qa_partial_36B7E875A19554DB711364BE0B072DA5.pdf)。RFF・run22090006–22121018に偏った部分統計。選択後Vzの範囲外/flow0、mixing Vz clamp0。4He候補4227（低R4108、高R＋TOF m²119）。PID再較正・純度・CF評価は行っていない。[Λ質量・アクセプタンス](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/preview_20261109/qa_page_03.png)／[4He β–η・TOF条件比較](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/preview_20261109/qa_page_07.png)。βは候補の|p|/Eであり、TOF実測βとは区別する。

14:33:56 JSTの独立Condor確認では再投入controller稼働、初期batch1–9の1800 jobs終了、metadata監査でROOT1575本、既知bad-run全入力/選択0が197 jobs、要確認28 jobs。28件はproducer完了footerなし、詳細原因は未確定。batch10は200 submitted、67終了、133 running、idle/held0。未投入5776入力（初期5725＋追加51）。終了メールは全batch後で未送信。再投入ROOTの全payload検証は後続hadd時に行う。

今回のquota実測：使用2.486998 TB / soft 2.588045 TB（96.10%）、hard 3.137801 TB、inDoubt 6.595 GB、inDoubt差引き空き94.452 GB。部分hadd終了でhadd予約0、再投入40＋安全20＋次200jobs5＝65 GBの投入条件を満たす。各batch前に再照会し、容量不足時はscheduler投入前にローカル待機する。holdは使用しない。

[確認記録JSON](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/auau7p7_2021_P22ib/status_20261109_001/status_check.json)。解析source/config/library・controller、schedulerは変更していない。ROOT/ログ削除、git commit/pushは行っていない。追加stderr診断の一部は自動承認review timeoutで未実行であり、28件の原因を確定したとは扱わない。


追記：保存した最新controller snapshotでは 2026-11-09T14:43:15.591813+09:00 に第10batchは終了82、running118、idle/held0。上記14:33の数字から進行している。quotaは 2026-11-09T14:44:15.404000+09:00 照会で、inDoubt 6.595 GBを差し引き、空き94.452 GB。必要65 GBを満たす。
