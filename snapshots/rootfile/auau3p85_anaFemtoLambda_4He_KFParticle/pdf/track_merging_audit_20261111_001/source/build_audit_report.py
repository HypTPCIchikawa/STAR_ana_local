#!/usr/bin/env python3
"""Assemble source, STAR-helix and saved-histogram evidence for the QA audit."""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);a=p.parse_args();out=a.study.resolve()
for name,key,value in [('source_audit.json','status','verified'),('helix_sign_validation.json','status','passed'),('hist_analysis_validation.json','passed',True),('plot_validation.json','status','passed')]:
 assert json.loads((out/name).read_text())[key]==value
tex=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{amsmath,booktabs,hyperref,pdfpages}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}
\begin{document}
\begin{center}{\Large $\Lambda$--He track-merging histogram audit}\\[4pt]
3.85 GeV fixed target, 2021 P24iy \quad 11 November 2026 (JST)\end{center}
\textbf{Finding.} The absence of a peak at $(\Delta\varphi^*,\Delta\eta)=(0,0)$ is reproduced in the saved histograms.
However, their construction has confirmed problems. The existing maps cannot establish that close-track effects are absent.

\textbf{1. Wrong curvature sign, including the \texttt{EventField} variant.}
\texttt{FemtoLambdaLegacy.cxx:997--998} uses
\[\Delta\varphi^*_{\rm current}=\varphi_p-\varphi_{\rm He}
 +\arcsin(A_p)-\arcsin(A_{\rm He}),\qquad A_i=\frac{0.3q_iBr}{2p_{T,i}}.\]
A direct SL24y \texttt{StPhysicalHelixD} position-at-radius calculation instead gives
\[\varphi(r)=\varphi_0-\arcsin\!\left(\frac{0.299792458qBr}{2p_T}\right),\]
for an outgoing track starting at the origin; $B$ is in T, $r$ in m, and physical $p_T$ in GeV/$c$.
The required curvature terms in the difference therefore have the opposite signs.

A ROOT 5.34/38 test of 192 tracks with both field polarities and charge signs verified this:
\begin{center}\begin{tabular}{lc}\toprule
Formula versus STAR helix spatial azimuth & Maximum difference [rad]\\\midrule
Current positive asin & 0.8665705\\
Negative asin with conventional coefficient 0.3 & 0.0003201\\
Negative asin with coefficient 0.299792458 & $1.33\times10^{-15}$\\\bottomrule\end{tabular}\end{center}
An ideal proton--He pair placed at exactly the same azimuth at $r=1.4$ m has real
$\Delta\varphi^*\simeq0$, but the current expression gives $+0.17935$ rad.
This is not an overall sign convention for the axis: reversing an axis would leave zero at zero.
This test proves the formula issue; it does not regenerate the measured sample.

\textbf{2. Primary-fit requirement on the secondary proton daughter.}
Lines 987--988 require \texttt{p->pPt()>0} and use \texttt{pMom()} for the $\Lambda$ daughter.
SL24y \texttt{StPicoTrack.h:72,178} returns zero primary momentum when no primary fit is stored.
The same analysis correctly uses \texttt{gMom()} for daughter acceptance at lines 848--850.
The angular QA consequently uses a different subset and daughter momentum from the selected $\Lambda$ sample.

\textbf{3. Additional legacy-only charge/rigidity issue.}
The unsuffixed histogram uses a fixed $+0.5$ T and $Z=2$ with unscaled nuclear rigidity $p_T$.
For physical nuclear charge and momentum, $2/(2p_T^{\rm rigidity})=1/p_T^{\rm rigidity}$;
the old expression doubles the nuclear curvature term. \texttt{EventField} fixes this factor and uses the event field,
but retains issues 1 and 2.

\textbf{Further geometric limitation.} The radius is fixed at 1.4 m. The asin argument is clamped without checking whether
both tracks reach that radius inside the TPC. The FXT vertex near $z=2$ m and secondary-track origins make this relevant.
A valid close-pair QA should use global track geometry and a common reachable radius, including TPC endcaps.
\newpage
\textbf{Saved-data audit.} The actual 1CE runtime source, Maker and histogram configuration were compared byte for byte
with the current files; all five checked files match. These are production issues, not later local changes.
The merge is the verified partial 2021 P24iy sample, SUMS
\texttt{1CE33994716F7896423A80AB994BAB14}; 26,945/26,994 subjobs, 49 excluded.

\begin{center}\begin{tabular}{lrr}\toprule
$^4$He QA quantity & Legacy & EventField\\\midrule
All $\Lambda$ signal SE pairs, all $k^*$ & 2,641,987 & 2,641,987\\
Merging histogram Entries, including flow & 2,253,813 & 2,253,813\\
QA fill fraction of selected signal pairs & 85.3075\% & 85.3075\%\\
Counts in displayed square $[-0.25,0.25)^2$ & 58,100 & 57,813\\
Fraction in displayed square & 2.5779\% & 2.5651\%\\
Counts in origin square $[-0.005,0.005)^2$ & 19 & 18\\
Counts in origin square $[-0.01,0.01)^2$ & 89 & 80\\
Counts in origin square $[-0.02,0.02)^2$ & 368 & 320\\\bottomrule\end{tabular}\end{center}
Both axes have 200 bins over $[-0.25,0.25)$: width 0.0025.
The $^4$He map averages about 1.45 entries per displayed bin; its origin $\pm0.01$ square contains only 64 bins.
The EventField central-to-adjacent-square-ring density ratio is $1.00\pm0.13$ there, using count-only errors.
The narrow range naturally puts most unrelated pairs in flow; flow is retained, not discarded data.
The 85.3\% fill rate is QA admission, not evidence for 14.7\% detector pair loss.

All six histograms have total-bin sums equal to Entries. The two conventions have equal Entries for each species.
Physical $^3$He and $^3$He with $^4$He mass are exactly identical, bin by bin, as geometry requires.
Global maps, origin zooms and projections confirm the lack of a sharp excess in the current definitions.

\textbf{Scope.} These are raw SE signal-window maps, integrated over all $k^*$ and centralities, for
proton daughter versus nucleus only. There is no pion-daughter map, angular ME reference,
or $k^*$--centrality--rigidity angular information in these histograms.
Track merging can reduce close-pair reconstruction efficiency; a peak in raw SE is not required.
Splitting and merging are distinct effects; see the detector-effects discussion in
\href{https://arxiv.org/pdf/1507.06842}{ALICE, arXiv:1507.06842, Sections 2 and 3.3}.
No ALICE cut value is transferred to this STAR sample.

\textbf{Relation to the high-$k^*$ CF.} The physical pair histograms are filled before this QA function,
and the current configuration has \texttt{closePairEnabled:false}.
The wrong QA formula therefore does not directly cut existing CF pairs.
The actual effect of detector pair efficiency on that CF remains unmeasured here.

\textbf{Required next validation.} Recalculate the angular QA from global daughter/nuclear track geometry,
with the correct field sign and charge/rigidity convention, and test reachable TPC radii.
Include both $\Lambda$ daughters and matched SE/ME views with $k^*$ and centrality information.
Validate with STAR helix and a small unchanged-selection PicoDst sample before any full production rerun or close-pair cut.
Existing two-dimensional maps cannot be corrected by a single axis flip or fixed shift, and cannot recover omitted tracks.

Production source, cuts and ROOTs were not modified. The accompanying nine atlas pages show the existing, uncorrected
histograms: three pages each for $^4$He, $^3$He, and mapped $^3$He.
\includepdf[pages=-,landscape=true,pagebox=cropbox]{merging_histograms_20261111.pdf}
\end{document}
'''
(out/'source/report.tex').write_text(tex)
md='''# Λ–He track merging QA監査（2026-11-11 JST）

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
'''
(out/'README.md').write_text(md)
print('Wrote audit explanation, report.tex and README')
