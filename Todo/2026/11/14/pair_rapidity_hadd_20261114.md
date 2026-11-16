# EP1/EP8共通成功統計のhadd

[詳細・予定output・状態・configlog](../../../../analysisnote/20261114/pair_rapidity_hadd_20261114.md)。

- [x] 両系列4058 jobs全終端・キュー空とROOT/log確認。EP8未出力35、EP1未出力34を記録。
- [x] 共通4023 jobs/72424 Pico referencesへ揃え、EP1だけの3496は保持して非選択と記録。
- [x] 50 unit testsと実ROOT5 2本smoke（1835hist加算・定数・8TH3D）PASS。
- [x] tmux `he4_pairY10_hadd_20261114_001` 起動、両系列の入力検証開始を確認。
- [x] 全選択入力audit・hadd・加算検証・元ROOT整理を完了する。
- [x] 定数とTH3Dを検証し、EP1/8の微小SE差を記録して許容する（ユーザー方針）。厳密一致の要件を更新。
- [x] 最小内容の結果メールを1通送信し、local sendmail受付を確認。
- [ ] 共通未出力34groupとEP8の3496の原本/転送を調査し、回復後の全統計比較を検討する。現在は自動再投入しない。
- [ ] EP分割の有用性のCF比較後に、保留中81runの較正を検討する。

10:59 JST現在は実行中。全投入統計ではなく、共通成功sampleの部分統合である。


---

## Session update — 15:09 JST / 2026-11-14：hadd要確認

共通4023/4058 jobsの部分統計。状態 `attention_required`、メールはlocal sendmail受付済み。

- epd8: `None`（None）
- epd1: `None`（None）

詳細は `rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/hadd_state.json`、`sent.json`、各merge stateのresult.jsonを参照。


---

## 17:59 JST 要確認メール調査

[原因と次の調査](../../../../analysisnote/20261114/hadd_attention_diagnosis_20261114.md)。

- [x] 停止原因をTH1F加算丸めと特定。最終41partの4QAヒストをbin単位でROOT5再現、両系列PASS。
- [x] 元8046ROOT不変、staging保持、正式出力未公開/未archiveを確認。
- [x] 別のSE差を確認。footer差100group、うち70group/82Picoにinput TTree entries差、priority2121/1273/1144/1932。
- [ ] staged Picoのpayload/replicaと再構成再現性は参考課題として保留。微小差は今回のhadd完了条件から外す（ユーザー方針）。
- [x] native精度検証と許容済み微小SE差のsnapshot照合後、正式出力へ確定・元ROOT整理を完了。元の失敗記録と旧通知markerを保持。

今回、再hadd・本番再投入・メール再送は行っていない。81runは保留のまま。


## 18:41 JST ユーザー方針更新・復旧準備

約222万SE中9ペア程度の微小差は許容し、完全一致と入力payload原因調査は今回のhadd完了条件から外す。全56SE総量差を記録し、native精度/全誤差/追加ヒスト検証後、既存ROOTでhaddを確定する。入力原因追跡の過去未完項目は参考課題へ変更。新規本番jobは投入しない。


---

## Session update — 19:03 JST / 2026-11-14：hadd復旧完了

ユーザー方針に従い微小SE差を許容。元ROOTから中間・最終までの固定TH1F例外をnative加算で全bin再現し、他histの従来監査、全hist誤差の有限性、pairY定数/8TH3Dを確認した。共通4023groupsを両系列で正式出力へ確定し、元8046ROOTを既存archive helperで保存した。

- epd8: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_DA0DEB4347BD834F1096F8D8798F2133_merge.root`
- epd1: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_FD44A5EB1B0AA22C546F169E86482587_merge.root`

完了メールはlocal sendmail受付済み。新規解析jobやPico再処理無し。再現性の完全一致を主張せず、SE差と入力entries差の記録を保持。未較正81runは保留。詳細：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/recovery_finalize_20261114_001`。soft quota余裕は完了時210.803GB。
