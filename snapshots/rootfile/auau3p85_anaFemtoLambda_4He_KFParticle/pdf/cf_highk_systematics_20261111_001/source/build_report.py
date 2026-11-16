#!/usr/bin/env python3
"""Build an explanatory summary around the saved raw-CF atlas."""
import argparse, hashlib, json
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);a=p.parse_args();out=a.study.resolve()
r=json.loads((out/'results.json').read_text())
rows={x['series_id']:x for x in r['series']}
def row(family,selection,group='cent0_20',species='4He',mass='signal',stage='selected'):
 return rows['|'.join([family,species,stage,mass,selection,group])]
def val(x,window='high07_08',latex=False):
 m=x['metrics']['windows'][window];return ('%.4f %s %.4f'%(m['value'],'$\\pm$' if latex else '±',m['error'])) if m['value'] is not None else 'NA'
def line(label,x,y):return label+' & '+val(x,latex=True)+' & '+val(y,latex=True)+r' \\'+'\n'
selections=[('Inclusive','baseline','all'),(r'$|n\sigma|<0.6$','nsigma','abs06'),(r'$|n\sigma|<1.0$','nsigma','abs10'),(r'$R<1.0$','rigidity','Rlt10'),(r'$1.0\leq R<1.5$','rigidity','R10_15'),(r'$1.5\leq R<2.0$','rigidity','R15_20'),(r'$2.0\leq R<5.0$','rigidity','Rhigh'),(r'$2.4\leq m_T<2.6$','mt','mt24_26'),(r'$2.8\leq m_T<3.0$','mt','mt28_30')]
table=''.join(line(label,row(f,s),row(f,s,'cent20_80')) for label,f,s in selections)
controls=''.join(line(label,row('baseline','all',species=sp),row('baseline','all','cent20_80',species=sp)) for label,sp in [(r'$\Lambda$--$^4$He','4He'),(r'$\Lambda$--$^3$He','3He'),(r'$^3$He with $^4$He mass','3HeAs4He')])
mass=''.join(line(label,row('baseline','all',mass=m),row('baseline','all','cent20_80',mass=m)) for label,m in [('Signal low half','signalLow'),('Signal high half','signalHigh'),('Low-mass sideband','SBNeg'),('High-mass sideband','SBPos')])
tex=r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}
\usepackage{amsmath,booktabs,hyperref,pdfpages}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}
\begin{document}
\begin{center}{\Large $\Lambda$--$^4$He high-$k^*$ systematics}\\[4pt]
Au+Au fixed target, 3.85 GeV, 2021 P24iy \quad 11 November 2026 (JST)\end{center}
\textbf{Finding.} The high-$k^*$ decrease is reproduced by the new histogram production.
Tightening the nuclear TPC pull alone changes the central-sample ratio little.
The decrease depends strongly on nuclear rigidity and also on pair $m_T$.
This identifies useful conditioning variables; it does not establish the physical or detector origin.

\textbf{Definition.} All quoted values are ratios of counts integrated over $0.7\leq k^*<0.8$ GeV/$c$,
scaled by $M_N/S_N$ from $N=[0.40,0.50)$ GeV/$c$:
\[ C_W=\frac{S_W}{M_W}\frac{M_N}{S_N}. \]
They are raw $\Lambda$ signal-mass-gate CFs, without purity corrections.
The reference interval is a comparison convention; no physical plateau is established.
Each cut is applied separately to a saved two-dimensional histogram.
$R=|\mathbf p_{\rm nucleus}|/2$ is in GeV/$c$; $m_T$ is in GeV/$c^2$.

\begin{center}\begin{tabular}{lcc}\toprule
Selection & Centrality 0--20\% & Centrality 20--80\%\\\midrule
''' + table + r'''\bottomrule\end{tabular}\end{center}

For centrality 0--20\%, the $|n\sigma|<1$ minus inclusive difference is
$-0.00139\pm0.00222$ after accounting for nested-count overlap.
The $R<1$ minus inclusive difference is $-0.05026\pm0.00412$.
These are count-model errors, not calibrated event-level significance tests.
The inclusive 0--20\% ratio is $0.8651\pm0.0028$ in $0.8\leq k^*<1.0$ GeV/$c$.
The inclusive 0--20\% linear slope over $0.5\leq k^*<0.8$ is
$-0.3299\pm0.0148$ $(\mathrm{GeV}/c)^{-1}$, with the full normalization covariance.

\textbf{Interpretation limits.} The high-rigidity selection also requires TOF and a mass-squared window;
the low-rigidity selection does not require TOF. Their difference cannot isolate a TOF effect.
A tighter pull window can change the rigidity composition; the joint $R$--pull distribution is absent.
Rigidity-correlated contamination is therefore not excluded. No new physics cut is adopted here.
\newpage
\textbf{Species and mass controls.} The same high-$k^*$ decrease is weaker in the selected $^3$He control.
Changing its mass hypothesis does not reproduce the full $^4$He decrease.
The two $^3$He representations use the same candidates and are correlated.
\begin{center}\begin{tabular}{lcc}\toprule
Control & 0--20\% & 20--80\%\\\midrule
''' + controls + r'''\bottomrule\end{tabular}\end{center}
The low and high halves of the $\Lambda$ signal window show similar decreases.
The sidebands are closer to unity; their different kinematics and composition prevent interpreting them as an independent purity correction.
\begin{center}\begin{tabular}{lcc}\toprule
$\Lambda$ mass region & 0--20\% & 20--80\%\\\midrule
''' + mass + r'''\bottomrule\end{tabular}\end{center}
\textbf{Production reproducibility.} With the same TH2D projections, the old/new inclusive 0--20\% values are
$0.90781/0.90788$; for 20--80\% they are $0.89431/0.89421$.
The productions contain strongly overlapping events; their difference is not tested as independent statistics.
Changing the normalization range changes the vertical scale and does not remove the declining shape.

\textbf{Input and validation.} New SUMS ID \texttt{1CE33994716F7896423A80AB994BAB14}.
The verified merge contains 26,945 of 26,994 subjobs, with 49 excluded; it is a partial merge.
It contains 2,099,784,951 input and 2,066,742,256 selected events.
ROOT 5.34/38 read 1,656 $k^*$ TH1/TH2 histograms from each production without changing either input.
The scan covers 1,680 series, 5 grouped plus 9 fine centralities, and 1,176 projection/partition checks.
Independent checks reproduce 318 numerical results and 24 finite-difference gradient checks.
Legacy TH2F/TH2D differences are confined to the $^3$He mixed-event $k^*\geq1$ overflow;
the maximum is 31 counts, $2.52\times10^{-7}$ relative, outside the plotted analysis range.

\textbf{Uncertainties.} Stored count variances are propagated through count integration, normalization reuse,
and overlap of nested cuts. Slopes use the full bin covariance from normalization.
The saved histograms do not recover event reuse, shared mixed-pair correlations, or model/purity systematics.
Missing-count bins are shown as unavailable, rather than zero.

\textbf{What this sample supports.} Independent pull, rigidity and $m_T$ projections, mass regions,
centrality and species controls are available.
Candidate storage and joint 8D Sparse were disabled for capacity reasons.
Simultaneous rigidity--pull cuts, low-rigidity TOF-tagged CFs, and arbitrary event-level remixing cannot be reconstructed from this file.
Pre-PID histograms retain track quality, species rigidity and TOF requirements and use the existing physical mixing pool;
they are diagnostic views, not an unbiased alternate production.

\textbf{Next discriminating measurement.} With candidates/PID held fixed, compare narrower rigidity or $m_T$ composition and
controlled mixing in multiplicity, vertex, run/field and event plane. Such a remix requires retained candidate records or a separate PicoDst rerun.
No new production is submitted in this study.

\textbf{Atlas guide.} Following pages: species (1--2), mass regions (3--4), pull (5--7), rigidity (8--10),
$m_T$ (11--12), pre-PID (13--14), normalization (15--16), and old/new production (17--18).
All configurations, counts and statistical checks are retained beside this report.
\includepdf[pages=-,landscape=true,pagebox=cropbox]{cf_highk_systematics_20261111.pdf}
\end{document}
'''
(out/'source'/'report.tex').write_text(tex)
md='''# 新1CEによる Λ–⁴He 高k*系統性（2026-11-11 JST）

新試料でも高k*低下は再現した。nσ窓を狭める効果は中央試料で小さく、rigidity・mT依存が大きい。原因の特定、最終CF規格化、PIDの新規採用は行っていない。

[結果と図付き報告](lambda_he4_highk_systematics_20261111.pdf)／[18頁atlas](cf_highk_systematics_20261111.pdf)／[全条件の数値表](metrics.tsv)／[計算・設定](results.json)／[独立検証](independent_validation.json)。

値はk*=0.7–0.8 GeV/cの積分raw CF、0.4–0.5 GeV/cで比較用規格化。純度補正なし。誤差は保存countに基づく統計誤差。

|条件|0–20%|20–80%|
|---|---:|---:|
'''
labels=['全選択⁴He','nσ ±0.6','nσ ±1.0','R<1.0','1.0≤R<1.5','1.5≤R<2.0','2.0≤R<5.0','2.4≤mT<2.6','2.8≤mT<3.0']
for label,(_,f,s) in zip(labels,selections):md+='|'+label+'|'+val(row(f,s))+'|'+val(row(f,s,'cent20_80'))+'|\n'
md+='''
0–20%のnσ±1とinclusiveの差は−0.00139±0.00222（nested overlapを考慮）。R≥2はTOF+m²窓を要求するため、R依存をTOF単独の原因としない。Rに相関する混入も、この独立1D投影だけでは除外しない。

³He物理質量／⁴He質量仮説の0–20%値は0.95828／0.95777で、⁴Heの全低下を質量変更だけでは再現しない。上下Λ信号窓は同様に低下し、左右sidebandはより平坦。旧04BFと新1CEの0–20%値は0.90781／0.90788。新旧はイベントが強く重なるため独立統計の差として有意性を付けない。

## 入力・検証

2021 P24iy、SUMS `1CE33994716F7896423A80AB994BAB14`。mainconf `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`。runmeta/configlog/merge resultは同IDで結合する。26,945/26,994、49件除外の部分統合。入力2,099,784,951、選択2,066,742,256 events。

新旧1,656histをROOT5で抽出し入力不変を検証。1,680系列、84,000表示bin、1,176投影・mass半窓検査。独立に318数値と24導関数を照合。³He ME overflowのみTH2F丸め差最大31count（2.52e−7相対）、0≤k*<1の結果へ影響なし。通常binの検査許容差は維持し、flowだけFloat ULPで別検査した。

[主検証](validation.json)／[新抽出検証](extraction_validation.json)／[旧抽出検証](extraction_validation_04BF.json)／[Float独立監査](legacy_float_independent_audit.json)／[描画記録](plot.log)。

## 再現

数値計算はホストのシステムPython/NumPyを使用する。CVMFS版NumPyはホストでlibgfortran不足なので、PYTHONPATHを除いて既存システム環境を使う。

```bash
env -u PYTHONPATH OPENBLAS_NUM_THREADS=1 /usr/bin/python3 STUDY/source/analyze_highk.py --settings STUDY/settings.json
```

抽出コマンドは`extraction_command.json`と`extraction_command_04BF.json`、入力指紋は`extraction_provenance*.json`に保存。描画は同mainconfを指定した`script/singularity_root_batch.sh ... root4star -l -b -q 'STUDY/source/DrawHighK.C+("STUDY/points.tsv","NEW_ATLAS.pdf")'`。新しい出力名を指定する。設定値は`settings.json`にあり、production設定を変更しない。

## 解釈と残作業

規格化0.4–0.5は物理plateauとして認定しない。count由来の規格化共分散とnested cut相関を含め、event/shared ME pair相関・純度やmassモデル系統は含めない。傾きは0.5–0.8をfull covarianceで記述し、低統計の欠測binは採用せず実際のfit binをmetricsに記録する。

candidate cacheとjoint8D Sparseがないため、R×nσの同時条件、低R TOF別CF、multiplicity/Vz/run/磁場/EPの再混合は未実施。prePIDは元physical poolのQAであり、quality/R/TOFを残す。³He core→⁴He混入率の定量や最終純度補正、controlled remixは未完了。

次の切り分けは、R/mT分布を合わせた比較と、候補/PIDを固定したcontrolled remix。後者には保存候補または別途PicoDst再処理が必要。今回、追加production・メール・git commit/pushは行っていない。

atlas頁: 1–2 species、3–4 mass窓、5–7 nσ、8–10 R、11–12 mT、13–14 prePID、15–16 normalization、17–18 production。先頭2頁を加えた報告PDFではatlas頁に2を足す。
'''
(out/'README.md').write_text(md)
print('Created report.tex and README.md')
