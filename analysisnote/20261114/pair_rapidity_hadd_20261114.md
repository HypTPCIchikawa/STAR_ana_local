# Λ–⁴He EP1/EP8共通統計のhadd

2026-11-14T10:59:00.458553+09:00 JST時点。ユーザー依頼のtmux実行・結合後メール通知を設定し、10:54 JSTに開始した。

## 対象と完了判定

各4058 jobsは全terminal・queue0、外側exit0は全件だった。ただしROOTと解析完了ログではEP8が4023成功/35未出力、EP1が4024成功/34未出力。EP1のみ成功した3496を比較用結合から除外し、元ROOTは保持する。共通4023jobs＝72424 Pico referencesだけを両系列で同じ番号・読み取り順序のgroupとして結合する。全投入統計のhaddではない。較正対象388runの入力から選び、未較正81runは引き続き保留。

[選択・欠落・未統合の詳細](../../rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/common_selection.json)、[固定job ID](../../rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/common_success_job_ids.list)。共通欠落34groupの原本/転送とEP8のみ未出力3496の回復は後続Todoで、本処理は新jobを投入しない。

## 実行と結果の場所

tmux: `he4_pairY10_hadd_20261114_001`、host: `starsub03.sdcc.bnl.gov`。driver pid3083961、EP8 pid3083987、EP1 pid3083988。両系列を並行し、入力ROOT全key再読→加算fingerprint→100本ずつROOT5 hadd→出力再読・加算一致→atomic publish→検証済み元ROOTだけsubjobs/IDへ移す。元ROOTデータは削除しない。

- epd8 SUMS ID: `DA0DEB4347BD834F1096F8D8798F2133`
  mainconf: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_covered25.yaml`
  runmeta: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_covered388p25_DA0DEB4347BD834F1096F8D8798F2133.json`
  configlog: `job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_covered388p25_DA0DEB4347BD834F1096F8D8798F2133.txt`
  予定output: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_DA0DEB4347BD834F1096F8D8798F2133_merge.root`
  状態: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_covered388p25_DA0DEB4347BD834F1096F8D8798F2133/status.json`、結果: `result.json`。
- epd1 SUMS ID: `FD44A5EB1B0AA22C546F169E86482587`
  mainconf: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_covered25.yaml`
  runmeta: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_covered388p25_FD44A5EB1B0AA22C546F169E86482587.json`
  configlog: `job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_covered388p25_FD44A5EB1B0AA22C546F169E86482587.txt`
  予定output: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_FD44A5EB1B0AA22C546F169E86482587_merge.root`
  状態: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_covered388p25_FD44A5EB1B0AA22C546F169E86482587/status.json`、結果: `result.json`。

[driver現在状態](../../rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/hadd_state.json)、[実行設定・固定SHA](../../rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/hadd_config.json)、[tmux開始証跡](../../rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/launcher_start.json)。現在は入力検証から進行中で、hadd完了やメール送信済みとは区別する。

## 検証と運用上の修正

- 共通IDを厳密に指定する`--include-job-ids`を追加し、SHA/選択理由をprovenanceへ保存。正常だが非共通の3496はoutside_requested_selectionで保持。
- 提出名にcovered388p25のlabelがあり実ROOT stemと異なるため、archive helperは凍結XMLのfromScratchを使うよう修正。runmeta/XMLを書き換えていない。
- archive unit20＋selection/merge unit30＝50 tests PASS、限定diff-check PASS。
- 実ROOT2本のROOT5 smokeで1835 histの加算fingerprint、TParameter定数、8TH3D・flow/誤差の有限性を確認。ROOT5にGetNcellsがない点を修正して最終macroを再確認済み。
- 全統計後にpairY nBins10/min−1.5/max0.5、He4自身Y threshold−0.5をYAML期待値と照合する。設定TParameterが4023倍にならないことも確認。
- さらに共通標本のhKstarSE* fingerprintがEP1/EP8間で一致することを確認する。CFの物理比較やnormalization判定は別工程。

## 容量と通知

作成前soft余裕264.693GB、tmux開始直前再確認264.339GB。現input約82.64GBは既に使用量へ計上済み。追加hadd staging保持の保守モデル約80GBより余裕がある。

通知先は既存の`yudai.ichikawa.d3@tohoku.ac.jp`。両系列のmerge・元ROOT整理・追加検証完了後に最小内容メールを1通送る。途中失敗は要確認メール。Message-IDはopaqueで、job名/ID/path/件数を本文に含めない。新規driver stateのlock/delivery_attempt/sentを使い、旧job完了メールとは別に二重送信を防ぐ。現時点では未送信。


---

## Session update — 15:09 JST / 2026-11-14：hadd要確認

共通4023/4058 jobsの部分統計。状態 `attention_required`、メールはlocal sendmail受付済み。

- epd8: `None`（None）
- epd1: `None`（None）

詳細は `rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/hadd_state.json`、`sent.json`、各merge stateのresult.jsonを参照。


---

## Session update — 17:59 JST / 2026-11-14：要確認メールの原因調査

[詳細調査](hadd_attention_diagnosis_20261114.md)。直接原因はhEpdStatus等TH1Fの加算丸め差が固定許容値1e-7を超えたこと。ROOT5で最終41partの4QAヒストをbin単位のFloat加算で再現し、両系列で全bin一致。元8046ROOTは開始時fingerprintから不変、最終staging各9.28GB保持、未公開/未archive。

別件として結合前SEにHe4 signal9/約222万ペアの差。入力一覧/順序とruntimeは同じだが、footer差100group中70group/82Picoで解析前TTree entriesが異なり、net入力差−8と一致。入力payload/replica差と再現性は未解明。正式公開・再hadd・新job・メール再送は行わず、診断JSONと調査対象groupを保存した。新規mainconf/jobid/configlog無し。


---

## Session update — 19:03 JST / 2026-11-14：hadd復旧完了

ユーザー方針に従い微小SE差を許容。元ROOTから中間・最終までの固定TH1F例外をnative加算で全bin再現し、他histの従来監査、全hist誤差の有限性、pairY定数/8TH3Dを確認した。共通4023groupsを両系列で正式出力へ確定し、元8046ROOTを既存archive helperで保存した。

- epd8: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_DA0DEB4347BD834F1096F8D8798F2133_merge.root`
- epd1: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_FD44A5EB1B0AA22C546F169E86482587_merge.root`

完了メールはlocal sendmail受付済み。新規解析jobやPico再処理無し。再現性の完全一致を主張せず、SE差と入力entries差の記録を保持。未較正81runは保留。詳細：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/recovery_finalize_20261114_001`。soft quota余裕は完了時210.803GB。
