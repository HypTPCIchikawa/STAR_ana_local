# oura／xjiang2のME実装とΛ–⁴Heの比較

調査対象はユーザー指定の`/gpfs01/star/pwg/oura/star-analyzer/`と、公開アカウント情報でXialei Jiangと確認した`xjiang2`の共有解析領域。`/gpfs01`は`/gpfs/mnt/gpfs01`と同じ場所。外部source/config/jobをread-onlyで調べ、解析・ビルド・ROOT読取・job投入は行っていない。個々のxjiang2配下のcodeの著者・最終production採用までを本人アカウントの帰属から推定しない。

## oura：保存production設定と下流実行

前回の手元checkoutだけでは欠落していたφ設定が、ouraの提出済みconfiglogから回収できた。

- SUMS tree producer: `629F1D1728AD40F21A4CA2E99A59FF68`。
- mainconf: `config/mainconf/main_auau3p85fxt_anaFemtoPhiTree_prod.yaml`。
- event Vz198–202 cmを実cut、横方向中心(-0.5,-2.2) cm、Vr<2.5 cm。
- saved mixing: Vz10区分、cent9=9区分、EP1区分、buffer20、bufferAll、cap0、bothDirections=true。
- 9/19全統計`phiqa_full_20260919`と9/28`r07_momnsigma25_notofopen06_v2_20260928`のblock0000 jobは、下流wrapperへbuffer=-1（設定値使用）、mode=bufferAllを明示している。後者のvariationはkaon PID/opening条件で、mixing modeは同じ。
- 9/19と9/28の凍結下流sourceも読めた。現行sourceとはhashが異なるため、過去jobの説明は保存sourceを優先する。9/19凍結版のME両方向loop・unit fill・全accepted event保存・oldest popを確認した。
- φp/dは同じcent9×Vzの過去poolと両方向総当たり。全accepted event（候補空も含む）が枠へ入り、混合後に現在eventを保存する。固定event gap・満杯待ち・1/Npool補正なし。
- run/Bはkeyでなく、下流block/call間でpoolは共有しない。StFemtoPhiTreeMakerは候補/イベント保存だけで、MEは下流で作る。異event KKから作るphi_mix背景とCF分母も別処理。

設定の根拠：[提出snapshot](/gpfs01/star/pwg/oura/star-analyzer/job/run/configlog/config_auau3p85fxt_anaFemtoPhiTree_prod_629F1D1728AD40F21A4CA2E99A59FF68.txt:25)、[mixing設定](/gpfs01/star/pwg/oura/star-analyzer/job/run/configlog/config_auau3p85fxt_anaFemtoPhiTree_prod_629F1D1728AD40F21A4CA2E99A59FF68.txt:86)、[全統計run記録](/gpfs01/star/pwg/oura/star-analyzer/job/merge/phiqa_full_20260919/RUN.md)、[9/19実行argv](/gpfs01/star/pwg/oura/star-analyzer/job/merge/phiqa_full_20260919/job/downstream_phiqa_full_20260919_block0000.sh:5)、[9/28実行argv](/gpfs01/star/pwg/oura/star-analyzer/job/merge/r07_momnsigma25_notofopen06_v2_20260928/job/downstream_r07_momnsigma25_notofopen06_v2_20260928_block0000.sh:5)。

保存source：[ME loop](/gpfs01/star/pwg/oura/star-analyzer/job/merge/phiqa_full_20260919/runtime/tools/anaFemtoPhiTreeDownstreamV3.C:1665)、[全採用eventのpool更新](/gpfs01/star/pwg/oura/star-analyzer/job/merge/phiqa_full_20260919/runtime/tools/anaFemtoPhiTreeDownstreamV3.C:1765)。

制限：ouraの現在configとruntime/configはOS EACCESで読んでいない。20は提出済みproduction snapshotの確認値で、下流argvのbuffer=-1は設定依存。各下流runに独自のbuffer上書きが無いことまでruntime/config本文で照合できたとは扱わない。9/28は確認した系統性runの例で、全runの最新版認定ではない。

## xjiang2：Λ–軽い核に近いdΛ系

主対象は`/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/`。`dlambda_p24ia`・`dlambda_tmva`でも主要なpool仕様は同じ。

|項目|保存sourceの処理|
|---|---|
|centrality|countrefmultによる9区分、cent<2を除外|
|Vz|198≤Vz<202 cmを8区分、0.5 cm幅、範囲外reject|
|EP|保存値psi_1_EPD_4/2を8区分、無効/端binはreject|
|保持|最大50イベント|
|入替|満杯後、ランダムな1枠を置換。直近50に限定されず古いeventが残り得る|
|組合せ|現在のd/t/³He×保存された過去Λのみ、全候補積、逆方向なし|
|間隔・warmup|固定event gapなし。過去1件以上から混合、満杯待ちなし|
|重み|1pair1count、1/Npool・pair capなし|
|入場|解析treeのevent cutを通った行。Λゼロでも保存するが、上流treeの候補保有条件は未確認|

このsourceのrunXiは3=d、5=t、7=³Heを選び、現在のruntestはtの例。dirnameだけで最終dΛ本番の粒子/selectionとは断定しない。EP fieldの物理定義・較正品質や解析結果は今回検証していない。参考sourceを2021 P24iyのproductionテンプレートやdata入力として採用していない。

根拠：[pool定義](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.h:84)、[ME loop](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:874)、[ランダム置換](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:1239)、[Vz区分](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:1564)、[EP区分](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:1617)。

## xjiang2：φp保存code

`/gpfs01/star/pwg/xjiang2/jxl/pphi/2.Analysis/CorrAna.C`と、これをcompileするMakefile・sandboxへ送るJobSubmitter.xmlを確認。派生版も存在し、最終本番で使った版・結果との対応は未確認。

cent9×Vz6区分（-35～25 cm）で、実効kMix10をランダムslot置換する。現在eventのslot0と保存slot1～10を含むi<jの55 slot組を毎event混合し、buffer同士の組も繰り返しfillする。fixed gapなし、raw k*はunit count、ME開始条件はcentだけのcounter≥10でVz別でない。

独立に静的確認した注意点：currentをME前にslot0と保存slotの両方へ書き、i=0とその保存slotの組を除外していないので同一eventがMEへ入る経路がある。MEでevent identityを持たず裸track IDを比較するため、別eventの同番号pairも除外され得る。これは現存sourceに対する所見であり、研究者の最終解析全体の正誤・現在採用版を断定しない。衝突energyも未確定で、3.85 GeVのVz設定として流用しない。

根拠：[currentと保存slot](/gpfs01/star/pwg/xjiang2/jxl/pphi/2.Analysis/CorrAna.C:192)、[ME loopとID比較](/gpfs01/star/pwg/xjiang2/jxl/pphi/2.Analysis/CorrAna.C:442)、[Vz区分](/gpfs01/star/pwg/xjiang2/jxl/pphi/2.Analysis/CorrAna.C:555)、[Makefile](/gpfs01/star/pwg/xjiang2/jxl/pphi/2.Analysis/Makefile:12)。

## Λ–⁴Heと比較して次に切り分けること

- [ ] Vz198–202の採用条件と端pool占有。φとdΛ参照は実cut、Λ1CEは範囲外clamp。
- [ ] 同cent9内のmultiplicity/運動量分布と、候補保有eventだけをpoolへ入れるΛ1CEのevent重み。
- [ ] 固定FIFO20とrandom replacement50の違いを、深さと更新方法を別々に変えて比較する。実event age/gapも保存する。
- [ ] EP区分の有無を適切な角度・validity方針でcontrolled comparisonする。参照コードのEP8をそのまま移植しない。
- [ ] 単方向/両方向を分けて比較し、入力順・run安定性への感度を調べる。

どのコードにも固定event間隔を空ける条件は見つからなかった。今回の比較だけでΛの高k*低下原因を特定したとはしない。解析コードの変更や再投入は未実施、track-merging追加確認は保留を維持する。

[原資料の指紋と読取範囲](external_me_comparison_evidence.json)／[先のΛ1CE audit](lambda_he4_me_construction.md)／[手元checkoutだけの旧比較](phi_p_d_me_comparison.md)。
