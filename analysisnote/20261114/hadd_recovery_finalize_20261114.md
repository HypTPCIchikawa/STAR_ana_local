# EP1/EP8 hadd復旧・確定

**日付:** 2026-11-14 JST  
**対象:** 2021 P24iy / Λ–⁴He pairY10 / 共通4023 group  
**前提:** [要確認メールの調査](hadd_attention_diagnosis_20261114.md)、[元の実行記録](pair_rapidity_hadd_20261114.md)

## 完了結果 — 19:03 JST

両系列の共通4023 job分を検証し、正式な結合ROOTへ確定した。元8046 ROOTを`subjobs/<SUMS-ID>/`へ保存し、EP1のみの非選択3496は元の場所に保持した。完了メールはlocal sendmail受付済み。soft quota余裕210.803GB。

- EP8: `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_DA0DEB4347BD834F1096F8D8798F2133_merge.root`（9282162513 bytes）
- EP1: `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_FD44A5EB1B0AA22C546F169E86482587_merge.root`（9282201802 bytes）

native加算は各42段階・4064 input readsで全bin一致。全1835 histの従来監査・Entries厳密一致、全bin content/error/raw Sumw2、56SEの固定snapshot、両系列のpairY定数/8TH3DがPASS。He4 signalのnet差は−9、bin差の絶対値和33、最大bin差3。差は記録して許容し、Picoの再解析は行っていない。

実行はtmux `he4_pairY10_hadd_recovery_20261114_001`。初回は既定tcshによるredirect解釈でPython開始前に終了したため、attempt未作成を確認しbashを明示して18:46:43 JSTに起動。失敗したlauncher記録も保持した。処理完了は19:03:51 JST、通知・完了記録を含む最終確認もPASS。以下の準備時点の記述は履歴として保持する。

## 更新した方針

ユーザーが「次に必要な作業を実施」「約222万ペア中9ペア程度の小さな差は気にしなくてよい」と指示したため、微小SE差は数値を記録して許容し、hadd復旧を進める。入力payload差の原因追跡・完全一致は今回の完了条件から外す。未較正81 runは引き続き保留。

全56 hKstarSEヒストの総entries/sum差は0.01%未満。最大絶対差はHe3 fullの53/15144888、最大相対差はHe4左SB YLabBelowの1/11347（約0.00881%）。He4 signalは9/2221363（約0.000405%）。bin番号付きfingerprintを物理的な差の大きさとして扱わない。全binのL1差・最大差も確認して保存する。

## 復旧方法

新規作業フォルダ：`rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/recovery_finalize_20261114_001/`

1. 元ROOT、既存監査、stagingのfingerprintと固定設定を確認する。最初の失敗結果と要確認メール記録をbackupして保持する。
2. 加算監査の例外はEP8 hEpdStatus、EP1 hEpdStatus/hRefMultWeightの3件だけに固定。TH1F/class・軸・label・entriesを検証し、元4023ROOT→round0 41ROOT→round1の全段階を実際の順でnative TH1::Addにより再現し、全binの厳密一致を要求する。その他histは従来加算監査を維持、全histのEntriesは厳密一致。
3. 全1835histについてcontent/error/raw Sumw2の有限性・非負分散、軸の同一性を確認する。56SEのbin差を記録し、事前に確認したsnapshotと照合する。QAはstaging前後fingerprintとsource/report SHAで固定して再利用する。
4. pairYの10分割/−1.5〜0.5、He4 Y threshold−0.5、8TH3Dの形・値・誤差を確認する。
5. 元8046ROOTの不変を再確認後、既存最終stagingをno-clobber hardlinkで正式出力へ確定。既存archive helperで入力をsubjobs/SUMS-IDへ保存する。ROOTデータの削除はしない。
6. 両系列完了後、別の一度限定通知stateから既存宛先へ最小内容の完了メールを1通送る。旧要確認メールのsent記録は削除しない。

## 実装と検証

- `replay_float_rounds.C`：native加算の再現。ROOT5実機smokeで正常なMergeを通過し、bin改変をexit84で拒否した。
- `check_recovered_hists.C`：bounded-memoryで全histのbin誤差とSE差を検査。ROOT5のdirectory cycleを段階解決する。内部Sumw2も直接検査。
- `finalize_recovery.py`：backup→検証→no-clobber確定→既存archive→完了通知→記録。独立レビューを実施し、Entries厳密・SEsnapshot照合・raw Sumw2検査を反映。

この節の時点では準備・検証中であり、確定・archive・完了メールを実施済みとは扱わない。実行状態は`recovery_state.json`、固定仕様は`recovery_config.json`、開始証跡は`launcher_start.json`、完了時の出力は追記へ記録する。

**mainconf / jobid / configlog:** 新規解析投入なし。対象は元hadd記録のEP8 `DA0DEB4347BD834F1096F8D8798F2133` とEP1 `FD44A5EB1B0AA22C546F169E86482587` のmainconf/configlog。Pico再読み取り・新規較正・本番source/config変更なし。


---

## Session update — 19:03 JST / 2026-11-14：hadd復旧完了

ユーザー方針に従い微小SE差を許容。元ROOTから中間・最終までの固定TH1F例外をnative加算で全bin再現し、他histの従来監査、全hist誤差の有限性、pairY定数/8TH3Dを確認した。共通4023groupsを両系列で正式出力へ確定し、元8046ROOTを既存archive helperで保存した。

- epd8: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_DA0DEB4347BD834F1096F8D8798F2133_merge.root`
- epd1: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_FD44A5EB1B0AA22C546F169E86482587_merge.root`

完了メールはlocal sendmail受付済み。新規解析jobやPico再処理無し。再現性の完全一致を主張せず、SE差と入力entries差の記録を保持。未較正81runは保留。詳細：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/recovery_finalize_20261114_001`。soft quota余裕は完了時210.803GB。
