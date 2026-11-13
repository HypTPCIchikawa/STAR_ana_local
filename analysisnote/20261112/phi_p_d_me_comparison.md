# φ–p／φ–dとΛ–⁴HeのME実装比較

現checkoutのsourceと残存設定を調査。`StPiMaker`は見つからず、実在する`StPhiMaker`について確認した。φp/dの最終production runと保存設定は特定できていないため、以下は現行実装の確認と、明示した残存設定の値を区別する。解析・job・ROOTの再実行はしない。

## 処理の担当

- StPhiMakerはφのK⁺K⁻不変質量背景MEを作る。φ–p/d CFの分母ではない。
- φp/dのPicoDst直接macroは共通StFemtoMakerを生成する。φpとφdの混合規則は同じで、候補選択とchannelが異なる。
- StFemtoPhiTreeMakerは全採用イベントと候補track情報・mixBinを保存する。producerはφ_rot/φ_mix/SE/ME pairを保存しない。
- tree方式のφp/d CF分母はtools/anaFemtoPhiTreeDownstreamV3.Cで、同event K⁺K⁻から再構成したφと別eventのp/dを組み合わせる。異event K⁺K⁻で作るphi_mix背景はさらに別の処理。

## φp/d CF分母（直接・tree下流）

同cent9×Vz×EPのpoolを使う。下流は通常、producerの保存mixBinを使い、recomputeMixBin optionで再計算する。φ(current)×p/d(buffer)とφ(buffer)×p/d(current)を設定に従って両方向に組み合わせる。同じeventUIDは除外。

φ側は候補空を含む全accepted eventをpoolへ追加する。直接StFemtoMakerのbuilderは候補ゼロでもspecies map keyを作り、全ベクトル空を除外する条件はLambdaLegacyだけなので同様である。ME後にcurrentを保存し、保持数超過で最古を捨てる。固定event gap・最低時間差・満杯までのwarmup待ちはなく、1pair=1count、1/Npool補正なし。run/Bはpool keyに含めない。下流呼出/block/slice開始でpoolは空、内部のrun/file/producer job境界resetはない。

Vzはイベント選択でminVz≤Vz≤maxVzを実cutし、同じ設定範囲をpool区分へ使う。Λの198–202を区分定義だけに使い範囲外accepted eventを端へclampする方式と異なる。φの具体的Vz数値は参照YAMLが欠落しているため未確定。

`bufferAll`は全eligible pairを処理する。`randomSample`でcapが正ならchannelごと・current eventごとに両方向候補積から一様非復元抽選する。mass/close-pair判定は抽選後なのでfill数はcap未満になり得る。sampling逆確率重みなし。

残存`config/cuts/mixing/mixing_auau3p85fxt_anaFemtoPhiTree_rndsamp.yaml`はnVz10/nCent9/nEP1/buffer20/randomSample/cap500/bothtrue。これは比較用設定の確認値であり、最終φproduction採用値とは未照合。下流wrapper・runnerのmode引数既定はbufferAllで、YAMLより優先する。正のbufferSize引数もYAMLをoverrideするため、mainconf名だけで実効mode/depthを断定しない。下流ROOTにはmixingModeのTNamedがある。

## StPhiMakerのK⁺K⁻質量背景ME

pool keyはVzのみで、cent9/EP/run/Bを使わない。eventのVz範囲外は実cutする。K⁺またはK⁻のあるeventをpoolに保存するが、混合のcurrent/buffer双方には両charge候補を要求する。current K⁺×buffer K⁻の単方向。pool eventとkaonを乱数で復元抽出し、試行数はmin(NK⁺current×NK⁻current,1000)。独自実装なので共通mixingMode/bothDirections/cap設定は使わない。1pair=1count、1/Npoolなし、満杯待ち・固定event gapなし。混合後にpoolへ保存しbufferSize超で最古を捨てる。

## Λとの比較と未確定事項

固定event間隔を空けない点・満杯待ちしない点・unit pair countはφとΛに共通。主要な実装差はφの全accepted event入場とVz実cut、Λの候補保有event入場と範囲外Vz clamp。φが別の固定event gapで問題を回避しているとは確認されない。これらの違いがΛの高k*低下の原因かは未検証。

φtree方式では候補選択を維持してbuffer深さ・混合modeをPicoDst再読なしで比較できる。ただし固定event-gap scanのoptionは現行にない。Λ1CEのhist-only出力は候補保存がなく、この点で異なる。

標準tree mainconfが参照するg5match mixingと3.85 event YAML、φ直接用mainconfは欠落していた。job/run/configlog/runmetaのPhi記録、job/merge、旧root-tree-study記録も現checkoutで見つからず、snapshotsは19GeV φ設定のみだった。別energyの設定を3.85の値へ代用しない。最終実行の数値はROOT metadata・実行argv・保存runtimeが揃ってから確定する。

## ソース

- [StPhiMaker混合](../../StMaker/StPhiMaker/StPhiMaker.cxx:700)、[event cut](../../StMaker/StPhiMaker/StPhiMaker.cxx:512)
- [tree producer event cut](../../StMaker/StFemtoPhiTreeMaker/StFemtoPhiTreeMaker.cxx:669)、[MixBins](../../StMaker/StFemtoPhiTreeMaker/StFemtoPhiTreeMaker.cxx:1081)、[イベント保存](../../StMaker/StFemtoPhiTreeMaker/StFemtoPhiTreeMaker.cxx:1433)、[producer metadata](../../StMaker/StFemtoPhiTreeMaker/StFemtoPhiTreeMaker.cxx:1572)
- [下流mixing](../../tools/anaFemtoPhiTreeDownstreamV3.C:1862)、[全event pool保存](../../tools/anaFemtoPhiTreeDownstreamV3.C:1963)、[override](../../tools/anaFemtoPhiTreeDownstreamV3.C:1101)
- [下流wrapper引数](../../script/singularity_run_anaFemtoPhiTreeDownstreamV3.sh:50)、[残存randomSample設定](../../config/cuts/mixing/mixing_auau3p85fxt_anaFemtoPhiTree_rndsamp.yaml:13)
- [前節のΛ1CE audit](lambda_he4_me_construction.md)


## oura側の実行記録による補足 — 22:13 JST / 2026-11-12

上記の「設定不明」は手元checkoutに限る。その後ユーザー指定oura repoで提出済み設定・全統計runと後続runのargv・凍結sourceを確認した。標準設定はVz198–202実cut・10/9/1・buffer20・bufferAll/cap0/bothtrue。下流runはbuffer=-1（設定依存）とbufferAllを明示。[更新比較報告](external_me_comparison_oura_xjiang2.md)を参照。旧本文を履歴として保持する。
