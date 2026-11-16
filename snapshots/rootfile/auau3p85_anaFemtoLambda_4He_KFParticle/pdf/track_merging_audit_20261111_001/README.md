# Λ–He track merging QA監査（2026-11-11 JST）

**結論：既存ヒストの作成に問題がある。原点ピークがないことだけで、近接trackの影響がないとは判断できない。** 元productionは変更せず、ソース監査・STAR helix試験・保存ヒストの数値と図を作成した。

[説明＋9頁の比較図](track_merging_audit_20261111.pdf)／[比較図のみ](merging_histograms_20261111.pdf)／[ソース・投入snapshot監査](source_audit.json)／[STAR helix試験](helix_sign_validation.json)／[ヒスト数値](hist_audit.json)。

## 確認した問題

1. **EventField版も曲率補正の符号が逆。** `StMaker/common/FemtoLambdaLegacy.cxx:997–998`はφp−φHeに+asin(Ap)−asin(AHe)を足す。STAR SL24yの位置角はφ(r)=φ0−asin(0.299792458 qBr/(2pT))であり、曲率項の符号が逆。±B、±電荷、複数pT/φ/pzの192ケースで検証した。現行+asinの最大誤差0.86657rad、正確な係数の−asinは1.33e−15rad。真に位置が一致する理想ペアを現行式はΔφ*=+0.17935radへ移すので、単なる軸の正負規約ではない。
2. **secondaryのΛ娘陽子にprimary情報を要求。** 同:987–988でpPt()>0・pMom()を使用。SL24y StPicoTrack.h:72/178によればprimary fitのないtrackはpPt=0。同じ解析の娘acceptanceは:848–850でgMomを使用している。QAの対象・角度が元CFの候補と一致しない。
3. **旧名はHeの電荷とrigidityの整合も悪い。** 同:992–993は固定B=+0.5T、raw rigidity pTへZ=2を掛けて曲率を過大にする。物理電荷2と物理pT=2×raw pTではZが相殺する。EventFieldはこの因子とイベント磁場を改善したが、上記1/2は残る。
4. **到達可能な半径を確認していない。** r=1.4m固定、asin引数をclampするだけ。FXT Vz≈200cm、secondary originを考慮して共通のTPC内到達半径・zを評価する必要がある。

実投入1CE runtimeのsource/Maker/hist YAMLなど5ファイルは現行とbyte一致。`source_audit.json`にarchive/member/hashを保存した。track IDを配列indexとして誤用する旧問題は修正済みで、今回の原因ではない。

## 実ヒスト

対象は`true/hDphiDeta_proton_{4He,3He,3HeAs4He}`と各`_EventField`、計6個。Λ信号窓のSEだけ、全k*・全centrality積分、proton娘のみ。ME・π娘・k*/R/centrality別角度軸はない。200×200bin、両軸±0.25、幅0.0025。

|⁴He|旧名|EventField|
|---|---:|---:|
|全SE信号ペア（k* flow含む）|2,641,987|2,641,987|
|merging Entries（flow含む）|2,253,813|2,253,813|
|QA fill率|85.3075%|85.3075%|
|表示範囲内|58,100|57,813|
|表示範囲内／Entries|2.5779%|2.5651%|
|原点±0.005|19|18|
|原点±0.01|89|80|
|原点±0.02|368|320|

表示bin平均約1.45count、原点±0.01は64binに80count。拡大・投影でも鋭い原点過剰は見えない。EventFieldの原点±0.01／隣接square-ring密度比は1.00±0.13（count誤差のみ）。表示範囲が狭いので約97.4%がflowに入ること自体は自然で、ROOT破損や削除ではない。fill率85.3%はQAのadmissionであり、検出器が14.7%のペアを失ったという意味ではない。

全bin積分=Entriesを6個すべて確認。旧名/EventFieldのEntriesが核種ごとに一致。³Heと³HeAs⁴Heは軸・内容・誤差・Entriesが完全一致。質量仮説だけではtrack幾何が変わらないという不変量を満たす。

## 解釈・CFとの関係

mergingは近接pairの再構成効率低下として現れ、raw SEの原点ピークを必ず作るものではない。splittingとは区別する。[一次資料：ALICE arXiv:1507.06842 §2/§3.3](https://arxiv.org/pdf/1507.06842)。この解析へ他実験のcut値を移していない。

物理CFは`FemtoLambdaLegacy.cxx:901–921`で先にfillされ、QAは926行から呼ぶ。現設定は`closePairEnabled:false`。今回のQA誤式が既存CFへ直接cutを掛けているわけではない。実際の検出器pair効率が高k*低下へ与える影響は未評価。

## 次の検証

- globalの娘・核track幾何、正しい曲率符号・電荷/rigidity・実磁場、secondary originと共通TPC内半径で再計算。
- protonだけでなくπ娘も、SE/ME双方、k*/centrality依存が調べられる形で保存。
- STAR helixでの試験と小統計PicoDstで選択を固定した比較後に、本番再生成・close-pair cutを検討。

今の2Dヒストだけを定数shift／軸反転しても修復できず、除外した娘も復元できない。候補cacheだけではmerging QAを再生できず、元trackへ戻る必要がある。今回は本番ソース・設定・ROOTの変更、追加SUMS投入、メール、Git commit/pushは行っていない。

## 再現・検証

抽出：`extraction_command.json` / `source/extract.sh`。解析：`python3 source/analyze_histograms.py --settings settings.json`（引数は`--help`も参照）。描画：`draw_command.json`。ROOT5/SL24yのmainconfは`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`。

helix runnerは`source/RunProbePhiStarHelix.C`、入力は`helix_synthetic_settings.tsv`と`helix_closepair_settings.tsv`、結果は同名`*_result.tsv`。位置角と接線の運動量角を区別して検査した。数学係数0.3の近似差と符号の誤りを分離している。

[ヒスト検証](hist_analysis_validation.json)／[描画検証](plot_validation.json)。説明2頁＋図9頁を最終PDFにまとめる。入力指紋・各成果物hashは`completed.json`へ保存する。
