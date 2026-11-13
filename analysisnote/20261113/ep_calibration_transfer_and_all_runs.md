# EP補正の共有と全run解析の可否

## 現時点の結論

2026-11-13 JST。未対応81runを含めてEP解析を準備することは可能。前回の「81runには追加較正が必要」は、現在のEP表に対応係数が無いという実装上の状態を述べたもので、81runすべてについて個別の係数を必ず新しく求めるしかない、と断定する意味ではない。既存係数や期間共通係数を共有する方法は、対象期間で妥当性を検証できれば選択肢になる。現時点では、その共有の適合性は未検証であり、有効化していない。

全runを残す方針に沿い、共有補正を独立標本で検証し、適合しなければ前期の期間共通／run別係数を準備する。補正表へのrun対応は明示し、既存runへのevent run IDの書換えや、未知runを無条件に通す処理は行わない。補正を得られない／EP無効なeventが全て使えるという意味でもない。

## PicoDst内にある較正と、解析用EP補正の違い

SL24y `StPicoEpdHit.h:22–28,101–106`はgain補正済nMIPと時刻依存DBのbad channel処理を説明する。一方、実P24iy sampleの保存StreamerInfoでは、StPicoEvent version8の82fieldにpsi／plane／Qx／QyなどのEP fieldは無い。StPicoEpdHitにはTObjectに加えてmId/mQTdata/mnMIPを保存している。

したがって「PicoDst信号の較正が未実施」なのではない。今回のrecenter（run・centrality別の平均Qを引く）とshift（recenter後角のsin/cos平均で角分布を補正する）は解析で選ぶhit/track/重み/centralityに対応した別段階であり、補正済EP角はPicoDstからそのまま取り出せない。

前期のnative Pico `recoveries/20261003/input_cache/st_physics_22123040_raw_5000023.picoDst.root` には4027entriesのEpdHit branchと非空payloadがある。標準StChain/StPicoDstMakerによる先頭100eventのsignal存在確認も実施中／結果はstudyの `early_epd_presence/` 参照。これはEP補正の流用妥当性や全81runの可否を判定する試験とは別である。

## Centralityでのrun上限拡張との相違

`StRoot/StRefMultCorr/Param.h:529`のfxtmult表は、共通のmultiplicity境界・Vz/trigger/luminosity係数と適用run範囲を保持する。`StRefMultCorr.cxx:865`はrun範囲に入る表のindexを選ぶ。既存修正はstop 22129999→22179999であり、境界／係数を変えず適用範囲を拡張している。ローカルPROVENANCEには公式適用との説明があるが、今回その外部一次資料は確認できなかったため、その表の転用妥当性を新たに実証したとは扱わない。

現在のEPは共通表のrun範囲でなく、TProfile2Dの(cent,runIndex)、TProfile3Dの(cent,shiftTerm,runIndex)ごとに異なる係数を読む。`FemtoEpdEventPlane.cxx:54–60,166,184,205,209`参照。run表に81runを追加するだけでは補正cellは生成されず、1000binsの空き枠にも係数が入るわけではない。run表をsortし直せば元の係数対応も壊れる。

元jhan Recenter.cxx348–356/460–467が平均Qを収集し、Shiftpar.cxx400–421がrun別recenter後の角度から20項sin/cosを収集する。共通recenterへ変更した場合は、その基準でshiftも再推定するのが整合的であり、旧run別shiftの単純平均で済むとは限らない。

## 適用対象と既存ファイル探索

未対応81runは22121033–22125011（day121–125）、既存388runは22158032–22179022（day158–179）。約1か月離れた前期への期間外適用になる。検出器応答やbad channelが同じとは未確認で、後期内の安定性だけで前期への適合性は保証できない。

jhan1/3p0とxjiang2/jxl/tree_p24iyの限定探索で、そのまま使える前期EP recenter/shift係数は確認できなかった。Runbyrun_P24iyの現run.hは前期81＋22124003の82runを持つが、隣接qa.rootは509run bins・416非空binsで現表と一致せず、元生成時のrun表を復元せずにrun対応を決められない。加えてそのsourceはtile中心とweight和で規格化前のraw Qを保存するrun-only QAで、現補正表と定義が違う。qaShift_p24iy.rootもpion/kaon/protonのPID平均shiftでEP用ではなかった。これらは採用しない。

## 全runを解析するための検証手順

1. 前期PicoからEvent/Track/EpdHitだけでraw normalized Qとevent/run/calibration-centを収集する。Lambda KF、核PID、pair/CF計算はこの収集には不要。
2. 既存recenter＋shift一式を共有する候補を定義し、対象81runの独立標本でrun×cent別に残留平均Q・EPのFourier成分・8区分占有率・hit状態を確認する。全run合算が平坦でもrun別の歪みが相殺され得るため合算だけで合格としない。
3. 共有が適合しなければ、前期の期間共通またはrun別recenterを求め、その基準の角度でshiftを求める。統計不足runは、応答が安定していることを検証した期間単位の共有も検討できる。raw-Q cacheを保存すれば、stage2/独立QAのためにPicoを繰り返し読む必要はない。
4. 較正生成用と異なる標本で確認し、対応run・係数・SHA・定義を明示して全runのEP1/EP8比較へ進む。必要統計や処理時間のbenchmarkはまだ無い。

新較正の生成時には適用側とcutsを統一する。原jhan Recenter/ShiftparのTPC eta窓A(-2.4,-1.305)/B(-1.205,0)とXialei適用側のA(-2.4,-1.28)/B(-1.18,0)には差があった。現在は元適用側の組合せを移植しているが、新生成時に無検討で同じ不一致を引き継がない。

recenterと角分布flatteningの区別・最終角の等方性確認は[Poskanzer–Voloshin, Phys.Rev.C58 (1998), Sec.II](https://arxiv.org/html/nucl-ex/9805001v2#S2)にも説明がある。

## EP分布と8分割の確認図

既存出力にはhEpdPsi1VsCent（補正後物理Psi1×cent9）、hEpdMixingBinVsCent（Xialei式で計算したEP区分0–7×cent9）、hEpdStatus（valid／除外理由）がある。raw/recenter/shiftの全段階やrun別のQAは、この3本だけでは揃っていない。共有係数の検証用に追加する必要がある。

今回、run22159051の先頭10000 Pico event／valid9374から、[PNG](../../rootfile/femto_lambda_ep_transfer_review_20261113_001/epd_eight_bins.png)と[PDF](../../rootfile/femto_lambda_ep_transfer_review_20261113_001/epd_eight_bins.pdf)を作成した。左上は全centの補正後EP角と45度ごとのnominal境界線、右上は実際のEP区分の占有数、下段は元2D。全統計の図ではない。角分布の平坦性に合格したとする図でもない。

概念上はPsi1の−pi～piをpi/4幅で8等分する。元のFloat_t変換／端点規約は維持し、図の境界線はnominal表示とした。実EP binのhistはコードでの実分類による。hEpdStatusがあるから、全入力eventが角度図に入るとは限らない。

## 作業範囲

今回、既存較正を81runへ適用するsource/config変更・新較正生成・SUMS本番再投入は実施していない。調査・実Picoの小規模signal存在確認・既存ヒストの確認図作成・Todo更新を行った。以前のproduction implementationの47file SHAは最終確認する。


## 前期signal確認・調査完了の追記

未較正run22123040の先頭100eventをcompiled StMaker/StChainで確認し、100/100に有限で正のgain補正済East hitがあった。総保存hit13911、有限・正East hit11072、East row1–8かつnMIP≥0.3のhit3470、非有限East hit0。該当row/閾値のhitも100/100に存在。これはこの1runの入力存在確認であり、EP有効性・flow・係数共有・残る80runの検証ではない。[source/log/JSON](../../rootfile/femto_lambda_ep_transfer_review_20261113_001/early_epd_presence/result.json)。

前回実装47filesのSHAは全て不変。新EP図は投影積分9374・8classesの閉包、PDF1頁のGhostscript描画、PNG目視確認PASS。共有係数の適用試験と新係数生成は今後の作業であり、今回調査だけで共有を合格扱いしていない。


## ユーザー方針更新 — 18:56 JST / 2026-11-13

未較正81runの解析・共有補正試験・新較正作成は保留とする。まず較正済388runでEP分割の有用性を確認し、有用と確認できたら81runの較正を行う方針。前節の共有補正検証／追加較正案は再開待ちであり、現在の実行予定ではない。
