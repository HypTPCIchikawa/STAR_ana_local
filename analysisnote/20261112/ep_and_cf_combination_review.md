# EP定義と、centrality/EPを分割したCFの合成

## 結論と対象

分割して同じ区分内だけでMEを作った後、SE/ME countを合算して割っても分割の意味は消えない。足し合わせる操作は、異なる区分間の禁止したmixed pairを復活させない。ただし、区分ごとのME相対量が違う場合は、合成の重みと規格化を明示しなければinclusive CFの形を変え得る。CFの単純な50:50平均も、ペア数を反映したinclusive CFとは別の観測量。

今回の再計算は既存1CEのselected inclusive Λ–⁴He signal（背景差引き・純度補正前）、規格化H=[0.4,0.5) GeV/c、cent9=6/7/8。Hは比較用の規格化域で、物理的plateauの認定ではない。元ROOTやproducer設定を変えず、既存exportから方式差を計算した。

## なぜ最後に合算しても分割の意味が残るか

二つのクラスA/Bを例に、分割混合した分母はM_AA+M_BB。最初から区分せず混合した分母はM_AA+M_AB+M_BA+M_BB。前者を加算してもM_ABやM_BAは加わらない。cent9、Vz、EPのどの区分にも同じ説明が成り立つ。

XialeiのdΛ producerは同cent9×Vz×EPでMEを作り、EP添字の無いhQ_Rebin_mix[cent9]とh2_QvsCent_mixへ直接fillしている。SEもEPを合算したhistに保存し、cent9は保持する。したがってEPを区分して混合し、EP合算countを保存する実例である。最終広中央度のCF後処理macroは今回見つからず、0–10と10–20の最終合成方法は断定しない。

根拠：[ME fill](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:1015)、[SE fill](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:689)、[hist定義](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:1377)。

## 合成する観測量と重み

S_i(k), M_i(k)を区分iのSE/ME、Hを共通規格化域とする。

\[
a_i=\frac{S_i^H}{M_i^H},\quad D_i(k)=a_i M_i(k),\quad C_i(k)=\frac{S_i(k)}{D_i(k)}.
\]

区分ごとに分母の規格化を合わせたcountから作るinclusiveな診断CFを

\[
 C_{\rm bal}(k)=\frac{\sum_i S_i(k)}{\sum_i D_i(k)}
\]

と定義すれば、代数的に厳密に

\[
 C_{\rm bal}(k)=\sum_i w_i(k) C_i(k),\qquad
 w_i(k)=\frac{D_i(k)}{\sum_j D_j(k)}
\]

となる。よって適切なk依存分母重みを使えば、区分CFの平均と、規格化を合わせたcountの和の比は一致する。

一方、規格化域のSEペア数による固定重みw_i=S_i^H/ΣS_j^HでC_iを平均する定義もある。これはk依存重みのinclusive CFと一般には異なり、どの平均を報告するかを選ぶ必要がある。逆分散重みも測定精度を優先する別の選択で、母集団ペア組成と同一とは限らない。centrality依存がある場合、等重み・固定ペア重み・k依存分母重みを混同しない。

実際の例としてALICE ΛKは磁場polarityごとにCFを分けてfitし、表示時に規格化域のnumerator pair数で重み付けして合成している。これは中央度を必ずその流儀で合成するという規則ではなく、固定ペア重み平均を用いた実例。

[ALICE, Phys. Rev. C103 055201 (2021), Sec.III F](https://arxiv.org/abs/2005.11124)。

## 生MEを足して全体で一度規格化する場合

現在の方法はa_all=ΣS_i^H/ΣM_i^Hを使う

\[
 C_{\rm global}(k)=\frac{\sum_i S_i(k)}{a_{\rm all}\sum_i M_i(k)}.
\]

a_iが全区分で同じならC_global=C_bal。a_iが違い、区分の相対的なM_i(k)の形も違えば、両者は形まで変わり得る。区分別規格化もHで相関が十分小さい/比較基準として共通に扱えるという仮定に依存し、ME形状そのものの不一致を直す操作ではない。

説明用の仮想count例：HではS1=S2=100, M1=1000,M2=2000。別のk区間でS1=50,S2=200,M1=500,M2=4000とする。各区分のaは0.1と0.05で、各C_iは両区間とも1。それでもglobalはHで1、別区間で250/(4500×200/3000)=5/6になる。balancedは両方1。この例は実データではない。

## 現在のmacroの実装

- 直近資料source/analyze_highk.pyはqa()140–152でcent9ごとSE/MEとerror2を加算し、rebin4後、group全体のHで一度規格化して比を作る。0–20はcent9 6/7/8、0–10は7/8、10–20は6。区分CFの平均ではない。
- SE/ME overlayのDrawSEME.Cは物理TH2を中央度範囲でProjectionXしてからgroup normを一度求める。
- common/macro/checkHistAnaFemtoLambda.C231–252もProjectionX→Rebin→ME規格化→Divide。TH1::Divideの表示誤差は規格化係数を固定として扱う。
- 直近studyは規格化countの導関数とk-bin間共分散も計算しているが、event/共有mixed pairの相関は含めない。
- 旧tools/auau3p85GeV/draw_LambdaNuclearId_relamom_sum_CFsub.C等も中央度count合算後の比で、旧norm域・background/purity処理は直近資料と異なる。

補足訂正：LambdaLegacy producer自体はSE/ME histogramを書き、CF除算を行わない。過去記録で触れたproducer共通設定のnormQ=0.5–1.0を、このΛ経路に保存CFがあるという意味で扱わない。今回の図の規格化は後処理の0.4–0.5。

根拠：[直近analysis source](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/source/analyze_highk.py:140)、[標準QA macro](../../common/macro/checkHistAnaFemtoLambda.C:231)、[producer Finish](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:963)、[Legacy Write](../../StMaker/common/FemtoLambdaLegacy.cxx:1002)。

## 既存データでの合成方式差

|区分|SE(H)|ME(H)|a_i=SE(H)/ME(H)|
|---|---:|---:|---:|
|0–5%|70954|2671573|0.02655888497|
|5–10%|23856|835605|0.02854937441|
|10–20%|27548|954253|0.02886865433|
|0–10%|94810|3507178|0.02703313034|
|0–20%|122358|4461431|0.02742572955|

ME/SE量は0–5/5–10/10–20で約37.65/35.03/34.64と一定ではない。

|0–20%の窓積分CF|現在global|0–10/10–20を別々に規格化したbalanced|3 finecentを別々に規格化したbalanced|
|---|---:|---:|---:|
|0.7–0.8 GeV/c|0.90787687|0.90710778|0.90641450|
|0.8–1.0 GeV/c|0.86510443|0.86457090|0.86394925|

これは各窓でcountを積分してから求めたCFであり、CF点の単純平均ではない。今回の中央度合成の差は約0.0005–0.0015で、高k*低下を解消しない。EP依存やVz内組成の効果を評価した結果ではない。

## 誤差

固定された正規化重みで平均する場合はVar(C)=Σw_i²Var(C_i)+2Σ_{i<j}w_iw_jCov(C_i,C_j)。正しく独立した区分なら後者は0だが、重み自体をdataから推定する場合はその揺らぎも加える。count ratioにも同じ注意があり、規格化域と各k*binの重なり、共通a_iによりbin間共分散が生じる。

実用上は保存subjob/block単位のjackknife/bootstrapで、各replicaの規格化・合成を最初から計算し直し、covariance matrixを保存する方針が分かりやすい。poolがjob内で完結する条件を使い、SEとMEを別々に独立再標本化しない。merged histだけから全event相関を復元できるとは扱わない。

## EPを導入する経路

一般のEPはQ_n=(Σw cos nφ, Σw sin nφ), Ψ_n=atan2(Qny,Qnx)/n。一次と二次では周期が2πとπで異なる。atan2の計算と区分化自体は小さい処理だが、受容のrecentring/flattening、subeventでの角度品質確認、候補による自己相関回避が必要。[Poskanzer–Voloshin, PRC58 1671 (1998)](https://arxiv.org/abs/nucl-ex/9805001)。

現行writerにはpT重みTPC Q2を保存し、candidate tracksを除外できる試験実装がある。Qx/Qy/角度/track数/validityを保存するがunflattenedの補助量であり、較正済みEPD角とは同じではない。Λ直解析のnEPを8へ変えるだけではguardで停止し、現実装の経路はwriterで保存→eventPlaneSource:storedでreplay。同じ有効EP event集合でEPなし/ありを比較する必要がある。

Xialei dΛは保存field psi_1_EPD_4を2で割って8区分するが、生成元のharmonic/較正定義は下流codeだけでは確定していない。Ψ1を2で割ることがΨ2の再構成になる訳ではなく、そのままTPCQ2へ置き換えない。

現在の1CEはhist-onlyでcandidate/EP情報を持たず、既存全統計k*histにEP区分を後付けできない。EPを追加するなら再処理時にcent9×EP別SE/MEとpool占有・有効partner数・event ageを保持し、分割数と混合深さを分けて比較する。各EP CFが同じ形と期待するinclusive解析でも、acceptanceや混合量の整合は確認する。実際に角度依存CFを測りたい場合は分割CFをそのまま保持し、合成で消さない。

根拠：[Q2保存](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:406)、[EP guard](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:237)、[invalid角度](../../StMaker/StFemtoMaker/StFemtoMaker.cxx:482)、[pilotと制限](../../StMaker/StFemtoMaker/README.md:293)、[Xialei field](</gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:287>)。

今回、producer/configの変更・新job投入・EPによる全統計再解析は行っていない。


## 数値再計算の保存・検証

既存exportから6ヒストを切り出す独立runnerを専用studyへ保存し、元export SHA/不変・通常200bin/表示50bin・既存global値・共通aでのinclusive方式一致・個別CF1でもraw合成が歪むtoy・共有入力を考慮した方式差のJacobian誤差を検証、全PASS。元ROOTは開いていない。

|窓積分CF|現在global|2区分balanced|2区分CFを規格化域SEペア数の固定重みで平均|
|---|---:|---:|---:|
|0.7–0.8|0.90787687|0.90710778|0.90763710|
|0.8–1.0|0.86510443|0.86457090|0.86496663|

固定重みは0–10%が約77.5%、10–20%が約22.5%。50:50平均ではない。表の固定平均は、各窓を先に積分した区分CFの平均で、点列の単純平均ではない。追加global再規格化なし。点列をinclusive D(k)でH内積分した平均が厳密1になるとは課しておらず、その値も別記録している。誤差は独立native countモデルでnorm/weight推定を伝播するが、event/共有ME pair相関は評価していない。

[study説明](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_centrality_combination_20261112_001/README.md)／[設定](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_centrality_combination_20261112_001/settings.json)／[計算source](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_centrality_combination_20261112_001/source/compare_combination.py)／[窓数値表](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_centrality_combination_20261112_001/window_results.tsv)／[表示点](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_centrality_combination_20261112_001/display_points.tsv)／[検証](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_centrality_combination_20261112_001/validation.json)。
