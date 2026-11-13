# 7.7 GeV 部分hadd・QAと再投入状況


## 状況確認 — 2026-11-09T14:44:15.404000+09:00

旧ID `36B7E875A19554DB711364BE0B072DA5` の部分hadd・出力検証・QA作成の完了を確認。回収2090/5085 jobsを統合し、2026-11-08 16:53 JSTにQAと通知を完了した。22 checkpoints、1786 hist、全payload readbackと段階加算、event/TTree件数の検証がPASS。入力168,680,371、選択後44,203,571、Λ候補13,729,465。再投入分は未統合。

[QA PDF（10ページ）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/auau7p7_qa_partial_36B7E875A19554DB711364BE0B072DA5.pdf)。RFF・run22090006–22121018に偏った部分統計。選択後Vzの範囲外/flow0、mixing Vz clamp0。4He候補4227（低R4108、高R＋TOF m²119）。PID再較正・純度・CF評価は行っていない。[Λ質量・アクセプタンス](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/preview_20261109/qa_page_03.png)／[4He β–η・TOF条件比較](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/pdf/qa_20261108_001/preview_20261109/qa_page_07.png)。βは候補の|p|/Eであり、TOF実測βとは区別する。

14:33:56 JSTの独立Condor確認では再投入controller稼働、初期batch1–9の1800 jobs終了、metadata監査でROOT1575本、既知bad-run全入力/選択0が197 jobs、要確認28 jobs。28件はproducer完了footerなし、詳細原因は未確定。batch10は200 submitted、67終了、133 running、idle/held0。未投入5776入力（初期5725＋追加51）。終了メールは全batch後で未送信。再投入ROOTの全payload検証は後続hadd時に行う。

今回のquota実測：使用2.486998 TB / soft 2.588045 TB（96.10%）、hard 3.137801 TB、inDoubt 6.595 GB、inDoubt差引き空き94.452 GB。部分hadd終了でhadd予約0、再投入40＋安全20＋次200jobs5＝65 GBの投入条件を満たす。各batch前に再照会し、容量不足時はscheduler投入前にローカル待機する。holdは使用しない。

[確認記録JSON](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/star_local/workspaces/auau7p7_2021_P22ib/status_20261109_001/status_check.json)。解析source/config/library・controller、schedulerは変更していない。ROOT/ログ削除、git commit/pushは行っていない。追加stderr診断の一部は自動承認review timeoutで未実行であり、28件の原因を確定したとは扱わない。


追記：保存した最新controller snapshotでは 2026-11-09T14:43:15.591813+09:00 に第10batchは終了82、running118、idle/held0。上記14:33の数字から進行している。quotaは 2026-11-09T14:44:15.404000+09:00 照会で、inDoubt 6.595 GBを差し引き、空き94.452 GB。必要65 GBを満たす。
