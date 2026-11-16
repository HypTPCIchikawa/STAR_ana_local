# 7.7 GeV collider本番job投入・容量保護

状態：2026-11-06 08:12:46 JST、tmux `auau7p7_production_20261106`で本番投入を開始した。SUMS IDは`36B7E875A19554DB711364BE0B072DA5`、最初の最大indexは5084（予定5085 subjobs、予算5086以内）。現在はSUMSのスクリプト生成中で、全件のCondor投入完了はまだ確認していない。ユーザーが「jobを投入、投入後に容量overにならないか確認」を指示した。解析設定は準備済みのまま、batch集約を最大5→10 files/jobに変更し、作業ノードscratch20 GBを明示する。

## 容量予算（GB=10^9 bytes）

| 項目 | 追加確保 |
|---|---:|
| 7.7 GeV ROOT（model約55 GBに余裕） | 72 GB |
| log/err/SUMS生成物/archive | 10 GB |
| その他 | 3 GB |
| 旧3.85 GeV hadd残作業 | 60 GB |
| 安全余裕 | 20 GB |
| 合計 | 165 GB |

投入直前08:12:47 JSTの実quota空きは173,560,823,808 bytes（173.56 GB、blockInDoubtを差引済み）、必要165 GBに対して8.56 GBの追加余裕がある。旧haddの既存中間使用分を予約60 GBから差し引かず、保守的に数えた。今後7.7 GeVのhaddを同時自動開始せず、終了後に別途容量を確認する。modelは絶対上限ではないため、実出力量とquotaを監視する。

## 確認済み

- [x] 30511 unique P22ib/2021入力と既知bad-replica LFN除外。mainconf/physics配置と15libraryは準備時と一致。
- [x] 解析設定の変更なし。analysis_infoの変更はfilesPerJob5→10だけ。template input_count guardは元からparameter化されており、生成後10となる。
- [x] min6/max10、最大予算5086 jobs。minはbest effortなので実際の最初の最大indexと全生成listで確認し、超過時は投入を止めて再評価する。
- [x] 最も大きい10 source files合計8,005,618,082 bytes。runtime/output分を含めたscratch request20 GB。
- [x] SUMS1.10.52のsimulateSubmissionで16入力→8+8の2jobs、request_disk=20000000 KBを確認。実farm dispatchなし。
- [x] quota guardの10 mock testsと、launcherのreport名なしID取得・予算超過・ID曖昧時の保護testsがPASS。実MTA/scheduler操作はtestで行わない。実condor_qが空queueで空stdoutを返す場合を追加修正し、実環境の読み取りcheck v2もPASS。失敗した初回checkを保持した。
- [x] tmuxで一回限定投入を開始し、今回ID/最初の最大indexと容量guard稼働を確認した。guard PID1269648、08:13:18 JSTのphaseはwatching、空き173.52 GB、対象queueはまだ0、fullSubmissionConfirmed=false。空queueだけでは監視を終了しない。
- [ ] 全投入完了後、runmeta ID/anaName/mainconf、連番coverageと実Condor submitted数を照合する。
- [ ] 全job終了時にROOT/producer logも確認し、完了または要確認メールを送る。

## 運用

作業領域：`star-analyzer_auau7p7_2021_P22ib`。mainconf：`config/mainconf/main_auau7p7_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`。
canonical XML：`job/joblist/joblist_auau7p7_anaFemtoLambda_4He_KFParticle.xml`。

旧haddは変更しない。容量guardは今回IDのOwner/UserLogを完全一致で限定し、余裕低下時にidleだけhold、回復時に自分のHoldReason分だけreleaseする。running・旧jobは操作しない。guard状態はGPFSquotaに依存しない/tmpへ保存する。旧haddの残予約とin-flight出力量も閾値に含める。

投入後の通知はreport名に依存せず、`Writting files for process HEX_index`の唯一IDと保存runmetaを照合する。宛先は既存のyudai.ichikawa.d3@tohoku.ac.jp、最小メール、一回送信。自動再投入・7.7GeVの自動haddなし。

投入wrapper PID1268554、submit.sh PID1268610。主状態とログは`/tmp/yichikawa/auau7p7_production_20261106_001/`に保存し、新workspaceの`tests/production_launch_20261106_001/launch_state.json`にも状態をmirrorする。全生成が終わりrunmeta・連番list・Condor submitted数が一致した時だけfull-submission markerを保存し、終了メール監視へ移る。生成が予算5086件を超えた場合やguardが停止した場合は、この新規SUMS投入プロセスだけを止めて要確認通知する。旧haddやrunning workerは止めない。

08:16:05 JSTの実queue確認：600件投入済み（running31、idle569）、全件request_disk=20000000 KB。生成済みsubjob list602本は各6入力だった。全件投入は継続中。08:15:19 JSTのguardはwatching、対象400件、hold0、空き168.64 GB、旧hadd残予約59.89 GB＋headroom20 GB＋running予約0.525 GBに対して余裕あり。通常pilot logは17,049 bytes/14,046 events（GPFS割当262,144 bytes）で、5085jobs×5–6小ファイルの割当約6.67–8.00 GBからlog/job10 GB予算を再確認した。異常な大量error logまでの上限保証ではない。

08:17:20 JSTには900件（running52/idle848）、hold0を確認。投入後の見かけの空き低下を再確認したところ、実使用量ベースの空き173.17 GBに対してblockInDoubt予約8.49 GBがあり、差引空き約164.68 GBだった。したがって初期の低下をROOT実出力増加とは扱わない。実ROOT回収0、生成中のjob/run通常ファイル割当は約0.76 GB。guardは引き続きinDoubtも差し引いた保守的な空きを使用する。

He4 PIDはGit標準tableの初期設定で、P22ib用purity/校正が未確定であることは準備時の説明を継承する。8D/candidate cacheなし、通常1786histと2ledgersを出力する。Git commit/push・データ削除なし。

[準備・検証記録](collider77_preparation_20261106.md)。実行・容量の詳細は新workspaceの`rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/production_launch_20261106_001/`に保存する。

小さな運用コード・設定・検証17ファイル（65,373 bytes）を[個人submodule内の保存記録](../../../../workspaces/auau7p7_2021_P22ib/production_20261106_001/README.md)へ保存した。ROOT/log/runtimeは含めず、commit/pushはしていない。

## 終了後の調査（11月8日）

全5085件の投入は11月6日08:35 JSTに完了。その後11月8日09:51 JSTに要確認通知が受理された。exit_0=2131／removed2954、回収ROOT2090、全queue0。容量guardがholdした2954件を再認識できずreleaseせず、未開始のままSUMSのhold2日上限で削除された。ここで記したhold/releaseの設計どおりには動作していなかった。正常終了や全統計完了とは扱わない。外側exit0だがROOTがない41件の内訳も含めて調査し、回収分はRFFのみ・run22090006–22121018に偏ることを確認した。

ユーザー指示により、検査済み2090本の部分hadd・10ページQAを11月8日15:17 JSTに別tmuxで開始した。旧guardは歴史的記録として保持し、修正版copyのmock検証のみ実施、解析の再投入はまだ行っていない。以後の状態・詳細は[11月8日のTodo](../08/collider77_hadd_qa_20261108.md)を参照する。
