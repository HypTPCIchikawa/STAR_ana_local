# Λ–⁴He：1CE productionのMixed Event作成方法

投入済み2021 P24iy、SUMS `1CE33994716F7896423A80AB994BAB14`を調査。mainconfは`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`、保存設定は`job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14.txt`。

保存runtime bundle `job/run/runmeta/runtime_bundle_auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14.tar.gz`（SHA256 `e27bd395d761e412533aff0cbaa9f667a40bc1ec5f41836dd50f87f8b4247151`）をread-onlyで確認。StFemtoMaker、FemtoLambdaLegacy、provider、KfEventSelection、CentralityHelper、MixingConfig、FemtoMixingSampler、FemtoCandidate、HistManager、run/ana macro、mixing YAMLは保存版と現行が一致するため、以下の現行行参照は本番実装に対応する。

## 実装で確認した事実

|項目|投入済み1CEの挙動|
|---|---|
|混合単位|イベントごとに保持した選別済み候補の四元運動量。回転・再samplingなし|
|中央度|cent9一致。0–20%は0–5/5–10/10–20%の各poolで作ったhistを後で加算|
|Vz|198–202 cmを10区分、通常幅0.4 cm。範囲外は端binへclamp|
|Event plane|nEventPlaneBins=1、このΛ経路はpsi2=-1でEP分割・整列なし|
|保持・使用数|bufferSize=20、maxMixEvents=20。直近の保存イベントから使用|
|組合せ|Λ(current)×⁴He(buffer) と Λ(buffer)×⁴He(current)の両方向の全候補積|
|event gap|固定offset・最低event数差・時間差の指定なし。直前の同一pool eventも対象|
|warmup|poolに1件あれば混合開始。20件になるまで待たず、SEの対応したwarmup除外もなし|
|重み|1 pair=1 count。pool数・有効partner数・候補数の1/N補正なし|
|境界|job内メモリpool。job開始時は空。今回の非replay入力経路ではfile/run境界resetなし|
|追加pool key|run/磁場/Vx/Vy/同cent9内のraw multiplicity/時刻は使わない|

`bufferAll`なのでpair cap・乱数抽選は使わない。未指定のmaxMixedPairsPerEventのデフォルト値を実効的なcapとは扱わない。

20枠は全input/accepted eventの直近20件という意味ではない。Λまたはd/t/³He/⁴Heのいずれかを持つeventを保存し、全候補空のeventを保存しない。CollectNucleiは全4核種のbest PIDをstoreへ入れるため、Λ/⁴Heが無くd/tだけのeventも枠を消費する。対象channelで実際に混合へ寄与する相手数と20枠を区別する。Λ候補はsignalだけに限定せず保存し、signal/sideband区分はpair fill時に同じ条件で適用する。

現在eventはSE/ME作成後にpoolへ追加し、21件になれば最古を捨てる。同一eventのMEをFillPairで禁止する。共有track判定はeventIndexとtrackIndexの両方を使い、別eventでtrack番号が同じだけの誤除外はない。SE/MEは同じk*計算・質量・HeのZ=2運動量補正を使う。

Vz198–202はeventの採用範囲ではない。Imp5互換選択ではfinite PV等を確認し、中央度側のbad-run・pileup・補正/validcent判定も別に通す。範囲外のaccepted eventがあれば端poolへ入り、実効的に先頭はVz<198.4、最後はVz≥201.6のaccepted eventを含む。実際の範囲外占有率は今回未集計。

## 今回の図から未確定な点と次の検証対象

高k*でscaled MEがSEより上という観測だけでは、以下のどれが原因か確定しない。相手を遠いeventへ変えるだけで改善するとも未検証。

- [ ] 既存ROOTの`mixingVzClampedEvents`と`hCentralityVsVz`等で端poolの実占有を確認。`hVz`は採用前にfillされるため単独で採用後Vzの証拠にしない。hCentralityVsVzも最終AcceptCentBinより前なので適用選択を確認する。
- [ ] 使用pairのΔinput-event・Δtimestamp、run/Bの一致、pool占有・有効Λ/⁴He partner数を記録し、event gapとbuffer深さを別々に変えるcontrolled remixを比較。
- [ ] pool立上がりと候補数積によるevent重み、同cent9/Vz内のmultiplicity・kinematics組成差をSE/MEで確認。

`hMixSamplerQA`はeligible/filled/forward/reverse等の総数を持つが、pairごとのevent gap・pool ageは保存しない。1CEはhist-onlyでcandidate保存なしのため、このhistだけから新しいevent-gap条件でMEを再作成できない。今回は実装調査のみで、設定変更・job投入・mail・元ROOT変更は行わない。track-mergingの追加Pico再読はユーザー指示による保留を維持する。

## ソース参照

- [投入済みmixing設定](../../job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14.txt:328)
- [処理順](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:392)、[mixing bin](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:2847)、[mixing loop](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:2876)、[pool更新](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:3085)
- [全核種候補の保存](../../StMaker/common/FemtoLambdaLegacy.cxx:744)、[event整合・pair fill](../../StMaker/common/FemtoLambdaLegacy.cxx:867)
- [event選択](../../StMaker/kfparticle/KfEventSelection.cxx:203)、[k*計算](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:2587)、[event/track identity](../../include/FemtoCandidate.h:113)
