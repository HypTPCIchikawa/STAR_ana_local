# Xialei方式のEPD一次EP・8分割の採用方針

ユーザーはXialei方式のEP分割をそのまま採用する意向。対象の計算・較正・分割定義を上流まで追い、移植可能な材料を確認した。採用対象はEPD一次角の計算・補正と8分割。dΛ側の粒子PID、centrality/Vz cut、buffer50/random replacement等をまとめて移植する決定ではない。

## 元の定義が確定した点

- 生成source：`/gpfs01/star/pwg/xjiang2/jxl/tree_p24iy/StRoot_v2/Shift/Shift.cxx`。
- EPD East（hit ID<0）のAB＝実row1–8を使う。RingGroup Aはzero-based0–3、Bは4–7。
- nMIP≥0.3、重みmin(nMIP,2)。PVからRandomPointOnTileへの方位角でQ1を作る。TileCenterを代用すると同じ計算ではない。
- Qをweight合計で正規化し、run×centralityでrecentring、20項のshift補正、[-π,π]へwrap。
- 保存値psi_1_EPD_4はpsi1shift_EPD[4]、つまり補正済み一次EP Ψ1。TPC Ψ2ではない。
- 下流はfloatへ読んだpsi_1_EPD_4を2で割り、getPsibinの(-π/2,π/2)8区分へ入れる。元保存角では概念的に45°幅の8区分に相当する。境界はopen、範囲外/無効はeventを除外する。

理想的な実数演算ではbin=floor(8*(Ψ1+π)/(2π))である。実装で厳密に元codeへ合わせる際はFloat_tへの変換・Float_t bin幅・端点の比較も元関数と照合する。単に一次角を2で割っても二次EPの再構成にはならない。

[Q1/補正と保存](/gpfs01/star/pwg/xjiang2/jxl/tree_p24iy/StRoot_v2/Shift/Shift.cxx:718)、[RingGroup](/gpfs01/star/pwg/xjiang2/jxl/tree_p24iy/StRoot_v2/StEpdUtil/StEpdEpInfo.cxx:170)、[下流角度読込](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:287)、[分割関数](/gpfs01/star/pwg/xjiang2/jxl/correlation/dlambda/CorrAna.cxx:1617)。

## 較正の所在と確認状況

- `/star/u/jhan1/pwgdata/3p0/EP_P24iy/recenter/recenterpar.root`（3,827,598 bytes）
- `/star/u/jhan1/pwgdata/3p0/EP_P24iy/shiftpar/shiftpar.root`（40,788,318 bytes）

両方とも読み取り可能で、stream SHA256を保存した。ROOT object/key/軸・各runとcentralityの較正entryはまだ検証していない。sourceのrun.hには22158032等の2021 runが388個あり、重複なしを確認。run表が2021であることと、今回全入力runへの較正coverage一致は別の確認。

較正lookupのcentは、全isPrimary track数を数え、閾値[5,10,19,33,51,77,111,155,183]で作る。現解析の中央度をそのまま代入せず、この定義で較正用centを再現・照合する。sourceはcent不足を9へ入れるが、EP較正の有効binを無確認でfallbackしない。

元ShiftはVz198–202、Vr≤2、trigger/pileup条件に加え、TPC A/Bと全7 EPD群のQ1/Q2 nonzeroまで確認してからpsi_treeを書いている。AB角だけを計算して使う実装では元の有効event集合と異なる可能性がある。元EP cacheはrunId/eventIdでjoinされるが、現在macro指定test_v2.rootは見つからず、完成済み全統計cacheの所在は未確定。

## 現フレームワークへの必要な接続

- Pico readerでEpdHit branchを有効化し、要求時にbranch存在を検証する。現在LambdaRunSupportはEvent/Track/TrackCov/BTofだけを読む。
- 元のQ1・較正形式に合わせるEPD providerを接続する。STARの標準StEpdEpFinderを呼ぶだけで元Shiftと同一とは保証できない。geometry/PicoHit APIは再利用できる。
- EPD Ψ1の物理角、validity、source/harmonic、較正file指紋・run表・geometry/RNG情報を明示保存し、既存TPCQ2用eventPlane2へ偽装しない。
- 現GetMixingBinは非負π周期のpsi2を仮定するため、Xialeiの負のpsi/2を直接渡してはいけない。raw EPとpool用座標/明示binを区別する。
- 現EPguardはreplay＋storedのTPCQ2に限定される。source/schema/configを対応させ、無効EPの扱いを明示する。EPあり/なしは同じvalid-event集合で比較する。
- ROOT5/STAR batch-matched buildを使い、元sourceとの角度/bin対応、run/cent lookup、境界とinvalid扱いを検証してから統計を拡張する。

標準StEpdEpFinderには、較正file欠落時の継続やinvalid sentinelのPsi wrapがあるため、使う場合もfile/key/軸・Qのstrict事前検証が必要。補正なしへのsilent fallbackは採用しない。

RandomPointOnTileは乱数を使う。同じEP定義の再現と、元runの各eventの角度のbit単位一致は区別する。後者には元のEP cacheまたは同一乱数系列・呼出順・geometryが必要。

## 現時点の到達点

同じEP定義と8分割を採用できるsource・2021run表・較正fileの所在は揃った。残るのは較正内容/対象runの対応と有効event集合の照合、および現frameworkへの移植・pilot検証。今回は方式の確定と記録であり、producer/source/config変更やjob投入はまだ行っていない。

既存1CE hist-only出力にEPD hit情報は無く、全統計histへ後付け分割はできない。導入後のMEにはPico EPD hit再読または同じeventのEP cacheが必要。track-merging追加QAの保留は維持する。

[原資料指紋と確認範囲](xialei_epd_method_sources.json)。
