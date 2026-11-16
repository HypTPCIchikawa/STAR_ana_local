# Λ–⁴He pair rapidity10区分・81run保留・quota確認

## ユーザー方針

未較正81runの解析・共有補正試験・新較正生成は保留。まず既存較正の388runでEP分割の有用性を確認し、有用と確認できれば81runの較正へ進む。今回81runの処理は再開しない。以前の案は再開待ちとしてTodoへ追記済み。

## 追加したpair rapidity診断

物理LambdaとHe4の、CFに使う同じ保存P4を足して `YpairLab=(P4Lambda+P4He4).Rapidity()` を計算する。単体He4 YLabや2粒子のrapidityの単純平均、cached y、重心系へのboostとは別。Heの物理運動量にはZ=2が既に反映されている。

範囲の指定が無かったため、lab −1.5〜0.5を10等分（幅0.2）として明示YAMLに設定し、commentaryで仮定を伝えた。`pairRapiditySplit: {enabled:true,frame:lab,nBins:10,min:-1.5,max:0.5}`。数値はYAMLから変更できる。

SE/MEとも同じpairYの定義と区分を使い、full mass/signal/leftSB/rightSBの4領域を保存。物理He4だけで、He3とHe3AsHe4は含めない。既存He4単体YLab<-0.5/>-0.5の32ヒストは独立に維持する。新しいpair選択、mixing pool分割、CF規格化、close-pair cutは追加しない。

8本のTH3Dが正本：`hKstar{SE,ME}VsPairYVsCent_lambda_he4{,_signal,_leftSB,_rightSB}`。X=k*（既存と同じ200bins）、Y=pairY10区分、Z=cent9の9区分。有限の範囲外pairはY underflow/overflowへ保存、非有限pairYは追加診断だけ省く。inclusive物理SE/MEは維持する。10区分を各TH1として複製せず、容量を抑えながら任意の隣接区分を後から投影・合算できる。

ROOT metadataは`FemtoHe4PairRapiditySplitDefinition`、`he4PairRapidityLabNBins/Min/Max`。定数metadataと既存He4YLab閾値はROOT5 TParameterのmerge-firstで保存し、haddで10×job数や−0.5×job数になるのを防ぐ。counterの加算は維持する。HistManagerはtype:TH3Dを明示的に作るよう対応し、既存TH3F既定動作は維持した。

## 設定と後処理

通常の診断入口は `main_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_ylab_merging.yaml` と `..._diagnostics_pairY10_noEp.yaml`。

今後のバッチ用には、既存dirty設定を上書きせず、独立したanalysis_info/anaName/output/logを持つ `..._epd{1,8}_pairY10_covered25.yaml` を用意した。固定入力は[covered388_pico_files.list](../../rootfile/femto_lambda_pair_rapidity_validation_20261113_001/covered388_pico_files.list)、元2021P24iy masterの順序を保持して未対応runを除いた73054files/388runs、除外7928files/81runs。joblistはfilelist入力で生成しているのでcatalogのrun範囲だけに依存しない。共有係数や未知run fallbackは追加していない。

25files/jobで約2923jobs/系列を計画。EP1とEP8は同じ入力リスト・分割条件を使う。SUMS実生成後のworker filelistでも同じgroupingか照合する必要がある。XMLは `job/joblist/joblist_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_covered388p25.xml`。XML構造・maxFilesPerProcess=25・共通inputlist・別output先を検証する。**star-submit/Condorへの本番投入は行っていない**。

隣接mergeは元TH3を保全し、`FemtoLambdaQa::ProjectPairRapidity(*h3,"SE_y12",1,2,6,8)`でYbin1–2/cent9=6–8のraw k*を得る。normalYbinは1–10、flowは0/11。単なるY軸rebinで元metadataを書き換えるより、派生投影を別objectへ保存する。ROOT5のnativeProjectionX(0,0)は全Yを意味し、表示zoomを継承する場合もあるため、公開helperは全X軸を直接作りbin content/varianceを明示合算する。helperのUF・OF・可変X軸・表示range不変も検証した。

## 検証

- SL24y / sl73_x8664_gcc485 / ROOT5.34.38のbatch-matched build PASS。libStarAnaConfig/libStCommon/libStFemtoMakerを再構築。Legacy headerを参照するStFemtoMakerも再コンパイル。
- 回帰1603559チェックPASS（新pair rapidity関連1140649）。独立log((sumE+sumpz)/(sumE−sumpz))、10区分・境界・flow・両ME方向・質量域・cent・隣接merge・既存診断不変・13不正設定・TH3Dの2^24超精度を検証。実2worker ROOT mergeで定数metadata不変とcounts加算を確認。
- 公開projection helperの独立小回帰130チェックPASS。既知weightedfixtureでUF1/normal6/OF8/all15と分散1/20/64/85、Xflow、他cent/Zflow排除、implicitPoisson、可変Xedges、表示range不変を確認。
- 同じrun22159051の先頭10000Picoで正常終了、accepted9761・Lambda174。旧EP8の1707 pair-physics hist、40002230cells（flowを含むcontent/error/Entries/statistics）が完全一致。He4 signal SE13/ME403を保持。新8TH3のY分割＋flow閉包、cent別・隣接mergeの24240cellsもPASS。
- この標本の8view合計はnormal988、Y underflow0/overflow0。全mass/signalなどの重複viewを足した数で、独立なpair総数ではない。rangeを全統計へ適用した際もflowを確認する。
- 新QAは既存17頁＋pairYの16頁＝33頁。SE/ME×4質量域について10区分を6panel/頁で表示する。k*はQA表示だけrebin4、rawpaircount、追加図で背景差引きやCF推定はしない。小統計の空パネルもそのまま表示。

[study](../../rootfile/femto_lambda_pair_rapidity_validation_20261113_001/README.md)／[実Pico比較](../../rootfile/femto_lambda_pair_rapidity_validation_20261113_001/real_pico_closure.json)／[QA PDF](../../rootfile/femto_lambda_pair_rapidity_validation_20261113_001/pairY10_qa_validated.pdf)。初期のQA underflow表示はROOT5の0,0特殊挙動に由来し、実TH3のflow異常ではなかった。修正後は全8viewの表示も0。初回回帰edge±1ULPの仮定も検証側をROOTnative丸めへ合わせて修正した。

## Quota：現時点の空きと投入後の増分予測

mmlsquotaで`rcf.bnl.gov:gpfs:star-pwg`を直接照会。user yichikawaのsoft上限2.588TB、hard3.138TB、使用約2.302TB。inDoubt控除後soft余裕は286056251392B＝286.056GB（decimal）。file数633237でfile quotaは0/0＝制限無し。rhstar groupも制限なし、fileset全体には約214TBの余裕があり、userのblock quotaが制約。問い合わせ時点で本人Condor queueとhadd/ROOT/submit/watcherは0。

旧26945 hist-only worker ROOT合計160.245GB、merged11.019GBを固定hist overhead＋event成分へ分けて見積もる。現入力は388run/73054files。新ヒストの圧縮前array増分を最大3.983MB/workerとして加算し、各系列のhadd中間保持39.8GB、小ログ等の256KiB割当、ROOT丸め、misc2GBを含む。

|files/job|予定jobs/系列|EP8単独の増分|EP1＋EP8の増分|
|---:|---:|---:|---:|
|5|14611|約210.1GB|約420.1GB|
|10|7306|約130.9GB|約261.9GB|
|25|2923|約83.5GB|約166.9GB|

**25files/jobを容量面で推奨し、専用設定とXMLを準備した。** 2系列で約119.2GBのsoft余裕が残る。既存ROOT成分をさらに2倍にするstress caseでも約219.2GB、余裕約66.9GB。5files/jobで2系列を同時に残す計画は現在のsoft余裕を超える。

同一10000eventの実測は旧4039238B→新4057009B（+17771B）、GPFS割当は双方4194304Bで増分0。新8TH3のcompressed TKey合計13895B。これは低占有のpilotの実測なので、全本番出力をこの差だけで見積もらず上のraw-array/stress caseも使った。

25filesを全処理した実worker出力を測った数字ではなく、圧縮・占有率の非線形性がある予測。実生成job数、25filesのruntime/worker staging disk、投入直前の最新quotaは実投入時に確認する。旧Condorの15h条件はCPU比10%未満時のremove条件で、絶対walltime上限とは確認できない。Picoのstagingはworker SCRATCHで、user GPFSへ全Picoを複製する計画ではない。削除・quota変更・本番job投入なし。

[quota監査](../../rootfile/femto_lambda_pair_rapidity_validation_20261113_001/quota_review.md)／[推計model](../../rootfile/femto_lambda_pair_rapidity_validation_20261113_001/quota_projection_models.json)／[input manifest](../../rootfile/femto_lambda_pair_rapidity_validation_20261113_001/covered388_input_manifest.json)。


## 中断後の再開・最終確認完了 — 21:40 JST / 2026-11-13

ユーザー再開指示後、最終QA33頁のSE/ME全10bin表示を目視、独立source reviewで阻害要因なし、git diff --checkと実テストsourceSHA一致を確認。最終source/config/library/XML manifestおよびcompleted.jsonをstudyへ保存。最新quotaは21:36:51 JST時点のsoft余裕286.167GB、2系列166.901GB後119.266GB残、stress219.186GB後66.981GB残。81run保留、本番未投入を維持。
