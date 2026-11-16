# hadd要確認メールの調査

**日付:** 2026-11-14 JST  
**対象:** 2021 P24iy / Λ–⁴He / EP1・EP8 pairY10 / 共通4023 group  
**関連:** [元のhadd実行記録](pair_rapidity_hadd_20261114.md)

## 結論

15:09 JSTの要確認メールの直接原因は、hadd後の加算検証でTH1Fの丸め差が相対許容差1e-7を超えたこと。両条件とも最終段階haddまで終了しているが、公開前の検証で停止している。容量超過・ROOT読取エラーを示す失敗ではない。元ROOT8046本は開始時fingerprintから不変で、最終結合ROOTもstagingに残っている。

別問題として、同じgroup番号・Pico一覧・runtime bundleでも、結合前からEP1/EP8のSEが完全一致しない。stdoutの入力情報を調べ、少なくとも70 group・82 Picoで同名ファイルのTTree::GetEntries()の差があると確認した。EP分割後の選択や今回のhaddでは説明できない、入力段階の差である。payload/replica差や再現性の原因は未解明。

今回の作業は診断に限定した。元の失敗記録・通知済みmarkerを保持し、正式出力の公開、archive、新job投入、メール再送は行っていない。許容差を無条件に緩めて完了扱いにはしていない。

## 1. メールの直接原因：TH1F加算精度

両`result.json`のerrorは`Merged histogram is not the input sum: hEpdStatus;1`。

|条件・hist|入力bin内容の合計|hadd後の合計|差|相対差|
|---|---:|---:|---:|---:|
|EP8 hEpdStatus|1818989493|1818989236|−257|−1.413e−7|
|EP1 hEpdStatus|1818989509|1818989269|−240|−1.319e−7|
|EP1 hRefMultWeight|1818989509|1818989696|+187|+1.028e−7|

全1835 histogramの名前・schemaは一致。上記が許容差を超えたhistで、Entriesは各系列の入力合計と一致した。既存監査はEntries、bin内容の総和、bin番号付きの総和で比較しており、誤差二乗和を検証したものではない。

TH1Fはbin内容を単精度で保存するため、大きいbinを段階的に足すと丸めが入る。原検証は`math.isclose(rel_tol=1e-7, abs_tol=1e-7)`を型によらず使用しており、この規模の加算で期待される丸め差も検出した。

### ROOT5での独立再現

診断フォルダ：`rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/recovery_20261114_001/`

`probe_float_merge.C`を本番と同じSL24y/ROOT5.34/38のwrapperで実行。保存済みround0の41ファイルを実際の順序で読み、binごとにFloat_t加算を再現した。EP1/EP8それぞれでhEpdStatus、hRefMultWeight、hCentrality、hCentrality16の最終ROOTの全binが完全一致し、不一致0だった。

このprobeは**最終hadd段階の加算の再現**であり、8046本の元ROOTからround0をbin単位で再構築した試験ではない。入力総和からround0までの丸めもある。厳密に復旧する際は、その範囲も明示した精度検証が必要。

- `diagnosis_summary.json`：加算不一致とprobe結果。
- `epd{1,8}_float_probe.tsv`：各binのdouble合計、Float加算、実際の値、誤差情報。
- `epd{1,8}_float_probe.log`：ROOT実行ログ。
- `epd{1,8}_round0.list`：加算順序。

## 2. 元ROOTと出力の保全

- 開始時manifestに記録された元8046 ROOTのsize/mtimeNs/inode/deviceを独立照合し、変更0。
- 各条件の入力監査はDONE4023、最終staging監査はDONE1。
- 各42本のhaddログにROOT Error/Fatal/segfault行なし。
- EP8最終staging ROOT：9282162513 bytes。EP1：9282201802 bytes。
- 両系列とも正式`*_merge.root`は未公開、`subjobs/<jobid>`へのarchiveも未実施。
- 15:09の自動記録の`None`は、公開output/archiveがまだ無いことを指す。元データ消失を意味しない。

元のdriverには再実行保護があるため、そのまま起動し直さない。既存stagingと失敗記録を保持した復旧手順が必要。

## 3. 別件：同じ一覧でも入力entriesとSEが違う

全4058 groupの入力`.list`はbytes・順序とも同一。採用した4023 groupは同じ72424 Pico rows。runmetaのruntimeBundle SHA・サイズ、commit、ROOT/STARも一致している。

一方、結合前の`input_audit.tsv`でhKstarSE* 56/56にfingerprint差があり、50本はEntries自体に差がある。

|集計|EP1|EP8|EP8−EP1|
|---|---:|---:|---:|
|input/requested/completed events|1842342512|1842342504|−8|
|reconstructed events|1818989509|1818989493|−16|
|selected Lambda|35875222|35875236|+14|
|He4 signal SE pairs|2221363|2221354|−9|
|He3 signal SE pairs|10839499|10839520|+21|

He4のSE差は約4.05e−6（約0.000405%）で小さいが、同一イベントの厳密な比較を検証したことにはしない。

### ログによる差の絞り込み

共通4023 group×2系列のstdout末尾を全照合し、input/reco/Lambdaのfooter差がある100 groupに絞った。全てrequested=completed、status=completed/0。

その100 groupのうち70 group・82 Picoで、同じbasenameの解析前`[input-file] entries=`が違う。差のnetは−8、絶対差合計122、1ファイル最大4。`analysis/FemtoLambdaRunSupport.h:115–116`で、これはREADで開いたPicoDst TTreeのGetEntries()の直接記録と確認した。

例：group2121、`st_physics_22166035_raw_*.picoDst.root`

|file末尾|EP1 entries|EP8 entries|差|
|---|---:|---:|---:|
|6500007|26440|26437|−3|
|6500008|25318|25314|−4|
|7000005|28217|28216|−1|

この−8はEP/centrality選択より前の入力差。調べた100 groupのログにはread error/ERROR/TFile::Recoverは見つからない。同じURL・basenameでも同じpayloadだったことはまだ検証できていない。

残る30 groupはfile別entriesが同じでもreco/Lambda数が異なる。同entriesは同payloadの証明ではなく、payload/replica差と再構成の非決定性を未切り分け。StRefMultCorrのglobal RNG使用はsource上存在するが、入力entries差の原因としては扱わない。

100/70 groupのIDは診断フォルダ`se_input_difference_summary.json`へ保存。優先調査は2121、1273、1144、1932。1144/1932はinput/reco一致でLambdaが−3/+5であり、別の再現性調査に向く。

## 4. 次に必要な作業（今回未実施）

1. 差groupのstdoutが示すstaged Picoの出自・物理replica・入力bytes/metadataを追跡する。同名・同entriesだけで同入力としない。
2. input/reco/Lambda差100 groupと、SE各bin差があるgroupの関係を確認する。footer一致だけで完全一致としない。
3. 必要な最小範囲で固定した入力を用い再現性を試験する。原因未確定で全jobを自動再投入しない。
4. hadd復旧時は、TH1Fのnative加算精度を説明できる検証、各パラメータ/TH3D、入力不変、SE標本の条件を確認し、結果と限界を明示して公開/整理する。
5. 高統計QAのTH1D化や精度に対応した監査は将来修正案。現時点で本番source/configや監査許容値は変えていない。

**mainconf / jobid / configlog:** 新規解析投入なし。対象2系列の詳細は元hadd記録の通り。既存388 runの範囲を維持し、未較正81 runは保留のまま。


---

## Session update — 18:41 JST / 2026-11-14：微小SE差を許容してhadd復旧へ

ユーザーの指示により、約222万ペア中9ペア程度の微小差を許容する。全56SEの総量差は0.01%未満、入力payload追跡/完全一致は今回の完了条件から外す。[復旧計画と検証](hadd_recovery_finalize_20261114.md)に従い、既存最終stagingをnative Float全段階再現、全hist誤差QA、pairY定数/TH3D検証後に確定/元ROOT保存/完了メールへ進める。新解析job・81run較正は行わない。現在準備/検証中。


---

## Session update — 19:03 JST / 2026-11-14：hadd復旧完了

ユーザー方針に従い微小SE差を許容。元ROOTから中間・最終までの固定TH1F例外をnative加算で全bin再現し、他histの従来監査、全hist誤差の有限性、pairY定数/8TH3Dを確認した。共通4023groupsを両系列で正式出力へ確定し、元8046ROOTを既存archive helperで保存した。

- epd8: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd8_pairY10_DA0DEB4347BD834F1096F8D8798F2133_merge.root`
- epd1: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10/covered388p25/auau3p85_anaFemtoLambda_4He_KFParticle_epd1_pairY10_FD44A5EB1B0AA22C546F169E86482587_merge.root`

完了メールはlocal sendmail受付済み。新規解析jobやPico再処理無し。再現性の完全一致を主張せず、SE差と入力entries差の記録を保持。未較正81runは保留。詳細：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/femto_lambda_pair_rapidity_hadd_20261114_001/recovery_finalize_20261114_001`。soft quota余裕は完了時210.803GB。
