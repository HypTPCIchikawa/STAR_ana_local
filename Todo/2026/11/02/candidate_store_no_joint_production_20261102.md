# 候補保存あり・8D PIDヒストなしの本番投入

状態：少統計検証完了。本番SUMSのjob生成・投入処理中。

- [x] hNuclearPidPairJointだけを明示的に無効にする専用設定を作成。
- [x] SL24y batch-matched build。
- [x] 2021 P24iyの4,027イベントで全denseヒストと候補保存のclosure、Sparseなしを確認。
- [x] 2021 P24iy catalog用XMLを確認し、投入処理を1回だけ開始。
- [ ] SUMS全subjob生成・全件投入完了（永続workerが継続中）。
- [x] 全job終了・ROOT回収・解析完了ログの確認後に通知する永続処理を設定。
- [ ] 実際の全job終了・通知送信（将来の処理）。
- [x] 現在の生成ID・状態・将来のconfiglog/runmetaパス・日報・個人snapshotを記録。

候補保存Writeを有効にして、通常2D PID QAと候補・PID/KF・混合条件は維持する。8D定義は専用QA YAML内でコメント化し、jointSparse: falseを明示する。通知先はyudai.ichikawa.d3@tohoku.ac.jp。候補ROOTは元job単位のUUIDを保つため、今回は自動haddを行わない。

## 少統計実行

2021 P24iyの4,027イベントで解析正常終了、accepted3,858、Λ62、status=completed/0。全denseヒストと候補payloadのexact比較を実行中。生成XMLの実catalog queryはproduction_3p85GeV_fixedTarget_2021 / P24iy / daq_reco_picoDst / nFiles=all、filesPerJob=5。SUMSの2列filelistをpath-onlyへ正規化し、解析失敗exit codeを伝播するguardを確認した。

mainconf: `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWriteNoJointSparse.yaml`。

## 少統計の最終比較

全1,786 denseヒスト・59,146,854 cell（flow含む）の軸・title/label・内容・誤差・Sumw2・Entries・統計量が完全一致。候補3,858行、物理＋QA224,365候補・4,487,300 Float値とevent scalar/typed vectorもbitwise一致。4He信号SE9/ME251。新ROOTのTHnSparse数0、joint metadataなし、producer completed。元PID・KF・混合・dense QA・候補保存は維持した。

本番catalogは149,817 locationから重複LFNを除いて80,982 unique files、検証passed。SUMS完了までsubjob総数は確定扱いしない。初回tmux起動は既定shellのredirect差で実行前に終了し、submit intentなしを確認後、明示的Bashで起動した。実際のsubmit.sh呼出しは1回のみで、exclusive intentにより自動再投入を禁止している。

[比較検証](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_no_joint_20261102_001/smoke_validation.json)／[投入・監視の状態](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_no_joint_20261102_001/launch_state.json)／[専用設定](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWriteNoJointSparse.yaml)。

## 投入・通知の現在状態（2026-11-02T22:42:51.021494+09:00）

- 生成中のSUMS ID: `2AF731F2C75544890B4E858354DDE43D`。最終report/runmetaはstar-submit帰還後に確定・保存される。
- launch_state: `submitting`。全件投入済みとはしていない。SUMSの生成が約1件/秒で進むため、完了まで数時間かかる見込み。
- 実行host: starsub03.sdcc.bnl.gov、tmux session: `he4_store_nojoint_20261102_001`、launcher PID20954。
- 完了通知: yudai.ichikawa.d3@tohoku.ac.jp、poll300秒、全件終了後のsettle300秒。投入完了後に通知監視本体を自動起動する。
- 通知成功の判定: 全Condorjob終了・queue空・ROOTヘッダー回収・requested==completed/status0。全payload検査済みとはしない。
- 候補ROOTのUUIDを保つため自動mergeなし。自動再投入・通知重複送信なし。
- 実生成worker26993はNoJointSparse mainconfを参照し、実入力は2021/P24iyの3ファイル。2列list正規化guardあり。

将来保存される投入記録: `job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_2AF731F2C75544890B4E858354DDE43D.txt`、`job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_2AF731F2C75544890B4E858354DDE43D.json`。通知記録: `job/run/notify/auau3p85_anaFemtoLambda_4He_KFParticle_2AF731F2C75544890B4E858354DDE43D/`。現時点ではこれらの生成を完了扱いしない。

[永続処理の設定](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_no_joint_20261102_001/launch_config.json)／[現在状態（更新されるJSON）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_no_joint_20261102_001/launch_state.json)／[今回の検証・投入開始記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_no_joint_20261102_001/validation.json)。

## 実行中のsource固定

実生成cshはworker開始時に元projectのanalysis/config/lib/include/StMaker/StRootを参照し、scratchへ実体コピーする。SUMS sandboxは投入時の不変archiveではない。従来運用どおり、この6ディレクトリは全job終了まで変更しない。今回のsource/config/libraryは投入開始後に変更しておらず、launcherに記録した全fingerprintの一致を確認した。投入中のSUMS・共有symlinkの変更や再投入はしていない。
