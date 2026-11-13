# Λ–⁴He pair rapidity10区分・EP1/EP8本番投入

開始は2026-11-13 22:05 JST。ユーザーの本番投入・完了メール依頼により、2021 P24iy・既存EP較正済388runの73,054 Picoファイルを対象に2系列を投入中。未較正81runは引き続き保留。

|系列|mainconf|SUMS ID|
|---|---|---|
|EP8|`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_covered25.yaml`|`DA0DEB4347BD834F1096F8D8798F2133`|
|EP1|`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_covered25.yaml`|`FD44A5EB1B0AA22C546F169E86482587`|

EP1もEP8と同じEP-validイベントを要求し、mixingのEP区分数だけを変える。両方にHe4単体YLab分割、pairYLab10区分、修正済みtrack-merging診断を含む。[事前検証](pair_rapidity_implementation_and_quota.md)済みsource/libraryを固定。投入直前に同じコマンドでXMLを再生成し、旧XMLも保存した。

## 22:22 JST：投入進行中

22:20 JST時点でEP8は2,100 jobsが受理（100 running／2,000 idle）、EP1は2,200（143 running／2,057 idle）、hold無し。入力workerファイル生成と投入は継続中であり、全job投入完了とは区別する。生成済み共通366groupsは、入力と各group内の順序が完全一致。全groupの一致と73,054ファイルのcoverageは投入完了後にランチャーが監査する。

maxFilesPerProcess=25は上限で、SUMSが`minSize=18,maxSize=25`で分割し、実際のjob数は4,058/系列見込み（最大index4057）。事前の2,923/系列は単純ceilの見積もりであり、この実job数で上書きして評価する。[再計算JSON](../../rootfile/femto_lambda_pair_rapidity_production_20261113_001/quota_actual_worker_projection.json)では、2系列＋hadd保持込みの増分は保守191.490 GB、stress256.347 GB、開始時soft余裕285.998 GBからそれぞれ94.508／29.651 GB残。自動haddは起動しない。経験モデルであり、実ROOTサイズ・live quotaを回収後に再確認する。データ削除なし。

## 永続実行・完了通知

ホスト`starsub03.sdcc.bnl.gov`のtmux session `he4_pairY10_epcompare_20261113_001`内でランチャー実行。投入は各系列1回だけ。切断後も投入し、submit.shが保存するrunmeta／全入力coverageの検証後、既存`script/watch_job_notify.py`を各系列に1つ自動起動する。

通知先は`yudai.ichikawa.d3@tohoku.ac.jp`、`--email-content minimal`、poll/settle各300秒。各系列の全worker終端・キュー空・ROOT/log完了を確認して1通ずつ通知する。失敗がある場合は正常完了とせず要確認通知。二重送信防止のlockと送信状態を保持し、今回テストメールは送信していない。現時点ではsubmit進行中のため、watcherの実起動は未確認。

[ランチャー状態](../../rootfile/femto_lambda_pair_rapidity_production_20261113_001/launch_state.json)、[設定・固定SHA](../../rootfile/femto_lambda_pair_rapidity_production_20261113_001/launch_config.json)、各`epd8/submit.log`／`epd1/submit.log`にcommandと全stdout/stderrを保存。正式runmeta/configlog/joblistlogへのリンクは投入終了後に追記する。


## 22:33 JST：投入・完了メール監視の実起動を確認

両系列submit.sh return0、各4,058jobsの全スケジューリング完了。両系列とも18files×4,048jobs＋19files×10jobs＝73,054 unique filesで、元対象に欠落・重複・余分な入力なし。全4,058groupの入力とgroup内順序がEP1／EP8間で一致。[全group照合](../../rootfile/femto_lambda_pair_rapidity_production_20261113_001/grouping_comparison.json) PASS。global元リストとはSUMSによる順序変更があるが、比較系列間では揃っている。両ZIP718files内容同一、固定検証manifestの対象SHAも一致。

メールwatcherはEP8が22:30:22 JST（pid2646420）、EP1が22:30:39（pid2646507）に実起動し、各4058jobsを認識。`settings.json`／`monitor.lock`／初回`status.json`と生存processを確認済み。全処理終了時に各系列1通。まだ解析完了でもメール送信済みでもない。

- epd8: [runmeta](../../job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_covered388p25_DA0DEB4347BD834F1096F8D8798F2133.json)、[configlog](../../job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_covered388p25_DA0DEB4347BD834F1096F8D8798F2133.txt)。
- epd1: [runmeta](../../job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_covered388p25_FD44A5EB1B0AA22C546F169E86482587.json)、[configlog](../../job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_covered388p25_FD44A5EB1B0AA22C546F169E86482587.txt)。

[最終投入・通知設定JSON](../../rootfile/femto_lambda_pair_rapidity_production_20261113_001/submission_summary.json)を保存。投入後live user soft余裕は277.927 GB（開始285.998GB）。見積もりには既に生成した記録を含むので、開始時基準の予測と区別する。

### 初期入力読み取りエラーと残作業

確認したstderr274本（EP1 150／EP8 124）のうち、EP1の5workers（3885/3904/3908/3909/3914）でstaged Pico1本がZombie/no keys recoveredとなり、preflightで停止。EP8はこの確認範囲で0件。全件エラー監査ではない。[LFN別記録](../../rootfile/femto_lambda_pair_rapidity_production_20261113_001/initial_input_readability_errors.json)。

SUMSの100% copiedは存在確認だけで、xrdcp戻り値・ROOT内容の保証ではない。原本破損と転送失敗は未区別。内側解析が失敗しても外側wrapperはexit0となるため、Condor exit0だけでは成功としない。既存watcherがROOT出力・処理完了ログまで検査し、異常なら要確認メールを送る。これらworkerは回収後に原本/転送を調査し、同じgroupを用いて回復してからEP1/EP8比較する必要がある。現在の5件だけを理由に全jobを取消・再投入していない。
