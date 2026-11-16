#!/usr/bin/env python3
"""Build the illustrated Japanese report from frozen measured results."""
from pathlib import Path
import json
import math
import re
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=ROOT/'share/figure/20261115_lambda_he4_ep_rapidity_cf'
DATA=json.loads((HERE/'analysis_results.json').read_text())
CURVES=DATA['curves']


def result(selection='inclusive',dataset='EP8',cent='c0_20'):
    return CURVES['__'.join((dataset,cent,selection))]


def pm(value,error):
    if value is None or error is None:return '不可'
    digits=max(0,min(4,1-int(math.floor(math.log10(error))))) if error>0 else 4
    return '$%.*f\\pm%.*f$'%(digits,value,digits,error)


def tail(selection='inclusive',dataset='EP8',cent='c0_20'):
    x=result(selection,dataset,cent)['windows']['k080_100'];return pm(x['value'],x['error'])


def shape(selection='inclusive',dataset='EP8',cent='c0_20'):
    x=result(selection,dataset,cent)['shape_ratio'];return pm(x['value'],x['error'])


def interval(i,merged=False):
    first,last=(2*i-1,2*i) if merged else (i,i)
    low=-1.5+(first-1)*.2;high=-1.5+last*.2
    return '$[%.1f,%.1f)$'%(low,high)


PREAMBLE=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=15mm,top=17mm,bottom=18mm,headheight=12pt]{geometry}
\usepackage{fontspec}
\setmainfont{DejaVu Sans}
\newfontfamily\latinfont{DejaVu Sans}
\newfontfamily\japanfont[AutoFakeBold=2]{Droid Sans Fallback}
\XeTeXlinebreaklocale "ja"
\XeTeXlinebreakskip=0pt plus 1pt
\newXeTeXintercharclass\japanese
\newcount\jpchar
\jpchar="3000 \loop\XeTeXcharclass\jpchar=\japanese\advance\jpchar by1 \ifnum\jpchar<"3100\repeat
\jpchar="3400 \loop\XeTeXcharclass\jpchar=\japanese\advance\jpchar by1 \ifnum\jpchar<"A000\repeat
\jpchar="FF00 \loop\XeTeXcharclass\jpchar=\japanese\advance\jpchar by1 \ifnum\jpchar<"FFF0\repeat
\XeTeXinterchartoks 0 \japanese={\japanfont}
\XeTeXinterchartoks \japanese 0={\latinfont}
\XeTeXinterchartoks 4095 \japanese={\japanfont}
\XeTeXinterchartoks 4095 0={\latinfont}
\XeTeXinterchartokenstate=1
\usepackage{amsmath,amssymb,fix-cm,graphicx,xcolor,array,booktabs,fancyhdr,tikz}
\usetikzlibrary{arrows.meta,positioning}
\usepackage[unicode,colorlinks=true,linkcolor=blue!60!black,urlcolor=blue!60!black]{hyperref}
\hypersetup{pdftitle={EP分割・ラピディティ選択と高k CF},pdfauthor={STAR 3p85 解析資料}}
\definecolor{navy}{HTML}{18354B}
\definecolor{teal}{HTML}{007F83}
\definecolor{warm}{HTML}{FFF1E4}
\definecolor{pale}{HTML}{EDF5F7}
\color{navy}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\scriptsize STAR 2021 P24iy / $\Lambda$--$^{4}$He signal}
\fancyhead[R]{\scriptsize EP・rapidity比較 / 2026年11月15日 JST}
\fancyfoot[L]{\scriptsize raw CF / 規格化 $0.4\leq k^*<0.5$ GeV/$c$ / 誤差は統計近似}
\fancyfoot[R]{\scriptsize\thepage\ / \pageref{reportend}}
\renewcommand{\headrulewidth}{.3pt}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}\setlength{\emergencystretch}{2em}
\renewcommand{\arraystretch}{1.35}\setlength{\tabcolsep}{4pt}
\raggedright
\newcommand{\He}{\ensuremath{{}^{4}\mathrm{He}}}
\newcommand{\phead}[3]{\clearpage\pdfbookmark[0]{#2}{p#1}{\small\color{teal}#1\quad #3}\par\vspace{2mm}{\fontsize{21}{27}\selectfont\bfseries #2}\par\vspace{3mm}}
\newcommand{\lead}[1]{{\large\bfseries #1}\par\vspace{2mm}}
\newcommand{\sub}[1]{\par\vspace{2mm}{\bfseries\color{teal}#1}\par}
\newcommand{\figcap}[1]{{\fontsize{9.2}{12.5}\selectfont\color{navy!75}#1\par}}
\newcommand{\boxnote}[1]{\par\colorbox{pale}{\parbox{\dimexpr\linewidth-5mm\relax}{\vspace{2mm}\small #1\vspace{2mm}}}\par}
\newcommand{\warmnote}[1]{\par\colorbox{warm}{\parbox{\dimexpr\linewidth-5mm\relax}{\vspace{2mm}\small #1\vspace{2mm}}}\par}
\newcolumntype{P}[1]{>{\raggedright\arraybackslash}p{#1}}
\begin{document}
'''


def build_tex():
    chunks=[PREAMBLE];page_count=0
    plotdir='plots_cropped' if (OUT/'plots_cropped/c0_20__inclusive__high.pdf').exists() else 'plots'
    def page(title,kind='結果'):
        nonlocal page_count
        kind=kind.replace('%',r'\%')
        page_count+=1;chunks.append('\\phead{%02d}{%s}{%s}\n'%(page_count,title,kind))
    def text(x):chunks.append(x+'\n')
    def figure(name,width='0.94\\linewidth'):
        p=OUT/plotdir/(name+'.pdf')
        if not p.is_file():raise FileNotFoundError(p)
        # Imported ROOT PDFs omit Tr; isolate them from XeTeX AutoFakeBold's
        # fill-and-stroke text state left by a preceding Japanese heading.
        text('\\begin{center}{\\normalfont\\special{pdf:literal direct q 0 Tr}'
             '\\includegraphics[pagebox=cropbox,width='+width+']{'+plotdir+'/'+name+'.pdf}'
             '\\special{pdf:literal direct Q}}\\end{center}')
    def tables_start(columns):text('\\begin{center}\\small\\begin{tabular}{'+columns+'}\\toprule')
    def tables_end():text('\\bottomrule\\end{tabular}\\end{center}')

    page('EP・rapidityで高k*のCFは1に近づいたか','先に結論')
    text(r'\lead{EP8で改善した。ただし、0--20\%全体の問題は残る。}')
    text(r'$\Lambda$--\He のsignal窓raw CFを、同じ規格化域で比較しました。以下の$C_{\rm tail}$は$0.8\leq k^*<1.0$ GeV/$c$の窓内countを積分してから作った比です。')
    tables_start('P{5.8cm}cc')
    text(r'条件 & EP1：$C_{\rm tail}$ & EP8：$C_{\rm tail}$\\\midrule')
    for label,sel in [('0--20\\% 全体','inclusive'),(r'0--20\%：$y_{\rm He,Lab}<-0.5$','he4_below'),(r'0--20\%：$y_{\rm He,Lab}>-0.5$','he4_above'),(r'0--20\%：$-1.1\leq y_{\rm pair,Lab}<-0.7$','pair_merge02')]:
        text(label+' & '+tail(sel,'EP1')+' & '+tail(sel,'EP8')+r'\\')
    text(r'10--20\%：$y_{\rm He,Lab}<-0.5$ & '+tail('he4_below','EP1','c10_20')+' & '+tail('he4_below','EP8','c10_20')+r'\\')
    tables_end()
    text(r'\sub{全体では、1からの差と下降の両方が残る}')
    text(r'EP8の全体は$C_{\rm tail}=0.906$、1より約9.4\%低い値です。$0.8$--$1.0$と$0.5$--$0.7$の形状比は$T=0.936$で、高k*側の低下も続いています。')
    text(r'\sub{ラピディティによって、残り方が大きく違う}')
    text(r'He4の小さい$y_{\rm Lab}$側では差が小さくなります。ペアの$[-1.1,-0.7)$はほぼ平坦ですが、中心値は約2\%低く、普遍的な解決の証明とはしません。10--20\%・He4 Belowの限定標本は、現在の統計近似内で1・平坦と整合します。')
    text(r'\warmnote{測定ヒストの上端は1 GeV/$c$です。ここでの評価は\textbf{0.5--1 GeV/$c$の範囲}に限ります。1 GeV/$c$より高い領域や、漸近的な収束は判断していません。}')
    text(r'\figcap{誤差は規格化の共有相関を含むhistogram統計近似。event/pair再利用による相関・系統誤差は未評価。EP1/EP8間の差の有意度を独立標本として計算していません。}')

    page('何をそろえて比較したか','データ・定義')
    text(r'\lead{固定した規格化で、EPとY選択の影響を見る。}')
    text(r'2021 P24iyの較正対象388 runから、両系列に共通する成功4,023 groupを使用。入力一覧は72,424 Pico参照で、未出力35 groupは除外しました。未較正81 runは含めていません。微小なSE差は既に許容した方針を維持しています。')
    text(r'\sub{EP1とEP8}')
    text(r'同じEP有効条件のイベントで、pool keyをcent9 $\times$ VzとするEP1、さらにEast EPD一次EPの8区分を要求するEP8を比較。最大20 event・両方向ME等の設定は共通です。')
    text(r'\sub{ラピディティ選択の意味}')
    text(r'He4自身の$y_{\rm He,Lab}<-0.5$／$>-0.5$と、ペアの$y_{\rm pair,Lab}=y(P_\Lambda+P_{\rm He})$の10区分を別々に比較します。どちらも\textbf{SE/ME共通のペア選択}であり、Yで独立したME poolを作る処理ではありません。同一eventのペアが異なるY区分に入ることもあります。')
    text(r'\boxnote{pairYとHe4単体Yのjoint軸は保存していないため、この2条件を同時に課した20区分は復元していません。今回の図は、それぞれ独立した選択系列です。}')
    text(r'\sub{CFと規格化}')
    text(r'\[C(k^*)=\frac{S(k^*)}{M(k^*)}\frac{M_H}{S_H},\qquad H=[0.4,0.5)\ {\rm GeV}/c.\]')
    text(r'各条件で$S_H,M_H$を求めますが、領域$H$は全条件で共通。元200binを4binずつ足し、$\Delta k^*=0.02$ GeV/$c$としました。対象cent9のSE・MEをそれぞれ加算してからCFを計算します。')
    text(r'\figcap{$H$は従来との比較用に固定した領域で、無相関のplateauであると実証したものではありません。本資料の「1との差」は、この規格化での相対残留です。}')
    text(r'\figcap{主結果は0--20\%（cent9 6--8）。cross-checkは0--10\%（7--8）と10--20\%（6）。Λ signal窓1.10708--1.12484 GeV/$c^2$。背景差引・purity補正前のraw CFです。}')

    page('割り算前のSEとME：全体','0--20%')
    figure('c0_20__inclusive__spectra','0.94\\linewidth')
    text(r'\figcap{各パネルのSEとMEを同じ座標に重ねました。MEは$H=[0.4,0.5)$の積分がSEに等しくなるように高さをスケール。上段EP1、下段EP8。高k*側でSEがスケール後MEより小さくなる形が、CFの1未満に対応します。}')
    tables_start('P{5.0cm}rr')
    text(r'表示範囲内の数 & EP1 & EP8\\\midrule')
    for label,key in [('SE：$0\leq k^*<1$','totalSE'),('SE：規格化窓$H$','normSE'),('ME：規格化窓$H$','normME')]:
        text(label+' & '+format(result(dataset='EP1')[key],',.0f')+' & '+format(result(dataset='EP8')[key],',.0f')+r'\\')
    tables_end()

    page('全体のCF：改善しても下降が残る','0--20%')
    figure('c0_20__inclusive__cf')
    text(r'\figcap{青EP1、橙EP8。灰色帯が規格化域、破線がCF=1。規格化域の外側を評価します。}')
    figure('c0_20__inclusive__high')
    text(r'\figcap{高k*拡大。EP8でも$C_{\rm tail}='+r'0.906\pm0.003$、$T=0.936\pm0.002$。レベルが1未満であるだけでなく、範囲内の下降も残ります。}')

    page('He4自身のラピディティで差が変わる','0--20%')
    figure('c0_20__he4_below__high')
    text(r'\figcap{$y_{\rm He,Lab}<-0.5$。EP8：$C_{\rm tail}$='+tail('he4_below')+'、$T$='+shape('he4_below')+r'。全体より改善する一方、中心値は約3.2\%低い。}')
    figure('c0_20__he4_above__high')
    text(r'\figcap{$y_{\rm He,Lab}>-0.5$。EP8：$C_{\rm tail}$='+tail('he4_above')+'、$T$='+shape('he4_above')+r'。こちらでは1からの差と下降が大きく残ります。}')

    page('ペアラピディティ依存を一覧する','0--20% / 10区分')
    figure('c0_20__pair10__tail_summary')
    text(r'\figcap{縦軸は$C_{\rm tail}$。$y_{\rm pair}\simeq-1.2$の上方点は規格化窓のSEが10しかなく、収束判定には使いません。空の区分をCF=0や1で埋めていません。}')
    figure('c0_20__pair10__shape_summary')
    text(r'\figcap{$T=(S/M)_{[0.8,1)}/(S/M)_{[0.5,0.7)}$。$T=1$は2窓の積分$S/M$が等しく、$T<1$は高い側の窓で小さいことを表します。窓内の形や平坦性は図・GLSも併用して確認します。規格化定数は相殺されます。}')

    page('数値一覧：全体・He4二分・pairY10','0--20%')
    text(r'値は窓内countから計算したCF。$T$は高k*の形状比で、CF点の算術平均ではありません。')
    tables_start('P{4.1cm}cccc')
    text(r'条件 & $C_{\rm tail}$ EP1 & EP8 & $T$ EP1 & EP8\\\midrule')
    for label,sel in [('全体','inclusive'),('He4 Below','he4_below'),('He4 Above','he4_above')]+[(interval(i),'pair_bin%02d'%i) for i in range(1,11)]:
        text(label+' & '+tail(sel,'EP1')+' & '+tail(sel,'EP8')+' & '+shape(sel,'EP1')+' & '+shape(sel,'EP8')+r'\\')
    tables_end()
    text(r'\boxnote{pairYの$[-1.5,-1.3)$、$[0.1,0.3)$、$[0.3,0.5)$では規格化用SEが無くCFは不可。$[-1.3,-1.1)$は規格化用SE=10で低統計。1が誤差内に入ることだけを、回復の根拠にしません。}')
    text(r'\sub{読み取れる傾向}')
    text(r'$[-1.1,-0.9)$では1・平坦と整合しますが精度は約3\%。$[-0.9,-0.7)$は約2\%低く、傾きは小さい。一方、Yが大きい側では低下が強くなります。EP8はその多くを緩和しますが、全域を1へ戻す結果ではありません。')

    for start in range(1,11,2):
        page('pairYの高k* CF：区分%d・%d'%(start,start+1),'0--20% / 全10区分を表示')
        for i in (start,start+1):
            sel='pair_bin%02d'%i;r=result(sel)
            figure('c0_20__'+sel+'__high')
            text(r'\figcap{区分%d：'%i+interval(i)+r'。EP8の規格化SE='+format(r['normSE'],',.0f')+'。$C_{\rm tail}$='+tail(sel)+'、$T$='+shape(sel)+('。低統計につき精密な収束判定は不可。' if i==2 else '。')+'}')

    page('隣接区分をまとめる：ほぼ平坦な領域','0--20% / pairY [-1.1,-0.7)')
    figure('c0_20__pair_merge02__cf')
    text(r'\figcap{幅0.2の隣接2区分を、SE/MEの加算後に規格化。区分境界をdataに合わせて最適化したものではなく、事前に定めた隣接2binずつのmergeです。}')
    figure('c0_20__pair_merge02__high')
    text(r'\figcap{EP8：$C_{\rm tail}$='+tail('pair_merge02')+'、$T$='+shape('pair_merge02')+r'。ほぼ平坦ですが中心値は約2\%低く、この部分標本の結果を全体の解決とはしません。}')

    page('区分を広げても、Yの大きい側は下がる','0--20% / adjacent pairY merge')
    figure('c0_20__pair_merge03__high')
    text(r'\figcap{pairY $[-0.7,-0.3)$。EP8：$C_{\rm tail}$='+tail('pair_merge03')+'、$T$='+shape('pair_merge03')+'。}')
    figure('c0_20__pair_merge04__high')
    text(r'\figcap{pairY $[-0.3,0.1)$。EP8：$C_{\rm tail}$='+tail('pair_merge04')+'、$T$='+shape('pair_merge04')+r'。mergeで統計は増えますが、この形状差は消えません。}')

    page('中央度を狭めても、全体の残留は続く','centrality cross-check')
    for cent,label in [('c0_10','0--10\%'),('c10_20','10--20\%')]:
        figure(cent+'__inclusive__high')
        text(r'\figcap{'+label+'、inclusive。EP8：$C_{\rm tail}$='+tail(cent=cent)+'、$T$='+shape(cent=cent)+r'。いずれも1未満です。}')

    page('He4 Belowの整合性にも中央度依存がある','centrality cross-check')
    for cent,label in [('c0_10','0--10\%'),('c10_20','10--20\%')]:
        figure(cent+'__he4_below__high')
        text(r'\figcap{'+label+r'、$y_{\rm He,Lab}<-0.5$。EP8：$C_{\rm tail}$='+tail('he4_below',cent=cent)+'、$T$='+shape('he4_below',cent=cent)+'。}')
    text(r'\figcap{10--20\%の限定条件では1・平坦と整合。一方、より中心の標本には低下が残るため、普遍的な回復とは結論しません。}')

    for sel,label in [('he4_below',r'$y_{\rm He,Lab}<-0.5$'),('he4_above',r'$y_{\rm He,Lab}>-0.5$')]:
        page('He4ラピディティ選択後のSEとME','0--20% / '+('Below' if sel=='he4_below' else 'Above'))
        text(r'\lead{'+label+'}')
        figure('c0_20__'+sel+'__spectra','0.93\\linewidth')
        text(r'\figcap{MEは選択後の同じ規格化窓でスケール。選択はSE/ME共通です。He4単体Yで変えたのはペア標本であり、独立したY別ME poolではありません。}')
        tables_start('P{4.8cm}rr')
        text(r'表示範囲内の数 & EP1 & EP8\\\midrule')
        for name,key in [('SE：$0\leq k^*<1$','totalSE'),('SE：規格化窓','normSE'),('ME：規格化窓','normME')]:
            text(name+' & '+format(result(sel,'EP1')[key],',.0f')+' & '+format(result(sel,'EP8')[key],',.0f')+r'\\')
        tables_end()

    page('統計量と、分割の閉包を確認する','QA')
    text(r'\lead{空・低統計の区分を「CF=1」と扱わない。}')
    tables_start('P{3.8cm}rrr')
    text(r'pairY区分 & EP8 SE$_H$ & EP8 ME$_H$ & EP8 SE$_{[0.8,1)}$\\\midrule')
    for i in range(1,11):
        r=result('pair_bin%02d'%i);t=r['windows']['k080_100']
        text(interval(i)+' & '+format(r['normSE'],',.0f')+' & '+format(r['normME'],',.0f')+' & '+format(t['SE'],',.0f')+r'\\')
    tables_end()
    text(r'\sub{閉包は108件すべて合格}')
    text(r'通常k* binでは、pairY10区分の和、YのUF/OFを含む和、He4 Below+Above、隣接mergeの内容が元ヒストに完全一致しました。0--20\%＝0--10\%＋10--20\%も独立に確認しています。')
    text(r'今回、$0\leq k^*<1$のSE/MEでpairY範囲外のcountはありません。k* overflowのMEにはEP1で3、EP8で2 countの閉包差が残りますが、CFの対象範囲外です。この差の起源は今回直接検証していません。')
    text(r'\boxnote{pairYの有限範囲外は元出力のUF/OFへ保存されています。He4二分の境界$y=-0.5$は両側から除外されます。今回は対象範囲の閉包差0であり、これらがCF低下を作ったという結果ではありません。}')

    page('誤差と「1への収束」の評価方法','統計上の範囲')
    text(r'\lead{値が1に近いことと、形が平坦なことを分ける。}')
    text(r'\[C_{\rm tail}=\frac{S_{[0.8,1)}}{M_{[0.8,1)}}\frac{M_H}{S_H},\qquad T=\frac{S_{[0.8,1)}/M_{[0.8,1)}}{S_{[0.5,0.7)}/M_{[0.5,0.7)}}.\]')
    text(r'図と表の窓CFは、CF点を平均せずcountを先に積分して計算します。$T$では規格化が相殺されるため、定数の変更だけでは高k*の形状差を消せません。')
    text(r'\sub{規格化を共有するbin間の共分散を含める}')
    text(r'\[\frac{\partial C_i}{\partial S_j}=C_i\left(\frac{\delta_{ij}}{S_i}-\frac{I_H(j)}{S_H}\right),\qquad \frac{\partial C_i}{\partial M_j}=C_i\left(\frac{I_H(j)}{M_H}-\frac{\delta_{ij}}{M_i}\right).\]')
    text(r'このJacobianと保存分散から共分散を作り、一次傾向のGLSにも使用しました。下表の傾きは$0.5$--$1$ GeV/$c$の要約で、物理モデルのfitではありません。')
    tables_start('P{6.3cm}c')
    text(r'条件（EP8、特記なし0--20\%） & 傾き $(\mathrm{GeV}/c)^{-1}$\\\midrule')
    for label,sel,cent in [('全体','inclusive','c0_20'),('He4 Below','he4_below','c0_20'),('He4 Above','he4_above','c0_20'),('pairY $[-1.1,-0.7)$','pair_merge02','c0_20'),('10--20\% He4 Below','he4_below','c10_20')]:
        f=result(sel,cent=cent)['fits']['high'];text(label+' & '+pm(f['slope'],f['slope_error'])+r'\\')
    tables_end()
    text(r'\warmnote{誤差はSE/ME histogram countを独立とする近似です。event/pair再利用の相関、EP1/EP8の共有相関、系統誤差は未保存・未評価。EP1/EP8の差の有意度を独立誤差の二乗和で評価していません。低統計・empty binの良い$\chi^2$は収束の証拠に使いません。}')

    page('結論と次の検証','今回分かったこと')
    text(r'\lead{EP8は有用な改善だが、残留のY依存を調べる必要がある。}')
    text(r'\sub{1. 全体の問題は、まだ解決していない}')
    text(r'全体の高k* CFは約0.866から0.906へ上がり、下降も弱まりました。それでも1より約9.4\%低く、$0.8$--$1$の範囲にも下降が残ります。EP分割だけで解決した、とは結論しません。')
    text(r'\sub{2. 特定のラピディティでは、1に近くほぼ平坦になる}')
    text(r'He4 BelowやpairY $[-1.1,-0.7)$で残留は小さくなります。特に10--20\%・He4 Belowは現在の近似統計内で1と整合します。ただし、Yや中央度によって結果が変わり、選択した部分標本の回復を全体へ外挿できません。')
    text(r'\sub{3. 次は「残る領域」と「小さい領域」の条件差を追う}')
    text(r'今回のY分割は共通のペア選択です。Yごとに独立したME poolを作る効果はまだ検証していません。必要なら、その方式や運動学条件をさらにそろえる方式を明示的に比較することが次の候補です。Yとrigidity/mTのjoint情報の有無も確認してから、保存項目・再ペアリング範囲を決めます。')
    text(r'\boxnote{今回、新しい本番job、Pico再解析、ME生成コードの変更は行っていません。未較正81 runも保留のままです。規格化域を動かして1へ合わせる操作は行っていません。}')
    text(r'\sub{判断の範囲}')
    text(r'この資料はsignal raw CFの0.5--1 GeV/$c$診断です。高k*残留の起源をMEだけに帰属する結果ではありません。最終的な相関の物理解釈や系統誤差の確定、1 GeV/$c$超の漸近挙動の主張には進んでいません。')

    page('再現に必要なファイル','付録')
    text(r'\lead{元ROOTを保持し、後処理の条件と数値を保存。}')
    text(r'共通成功4,023 groupのEP1/EP8結合ROOTを使用。前日のhadd復旧結果、設定定数、入力不変性、元ROOT保存の記録を継承しています。')
    text(r'\sub{本番識別子}')
    text(r'EP1 SUMS ID：\texttt{FD44A5EB1B0AA22C546F169E86482587}\\EP8 SUMS ID：\texttt{DA0DEB4347BD834F1096F8D8798F2133}')
    text(r'\sub{後処理study}')
    text(r'{\small\texttt{rootfile/femto\_lambda\_ep\_rapidity\_cf\_20261115\_001/}}')
    text(r'\begin{itemize}\item \texttt{study.yaml}, \texttt{input\_manifest.json}：条件・入力fingerprint。\item \texttt{raw\_projections.root/tsv}, \texttt{source\_axes.tsv}：264 raw投影と元軸。\item \texttt{analysis\_results.json}, \texttt{metrics.csv}, \texttt{cf\_covariances.npz}：132 CF系列・定量・共分散。\item \texttt{cf\_histograms.root}：SE/ME/スケールMEと、有効binだけのCF graph。\item \texttt{closure\_checks.json}：108閉包。\end{itemize}')
    text(r'\sub{計算と表示}')
    text(r'ROOT 5.34/38・SL24yのbatch-matched wrapperで抽出/作図。数値はsystem Python3のisolated mode、NumPy 1.23.5を使用。規格化Jacobianは既存\texttt{analyze\_highk.py}の式を再利用し、別の式で独立再検算しました。')
    text(r'\sub{図とPDF}')
    text(r'{\small\texttt{share/figure/20261115\_lambda\_he4\_ep\_rapidity\_cf/}}')
    text(r'元の180 vector PDF、crop版、組版source、provenance manifest、ページ描画を保持。入力ROOTの原本は変更していません。')
    text(r'\figcap{この資料の数値・図は今回の実データから生成しました。理想CFや説明用toy dataを測定図へ混ぜていません。}')
    text(r'\label{reportend}\end{document}')
    return '\n'.join(chunks),page_count,plotdir


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tex,pages,plotdir=build_tex();(OUT/'cf_ep_rapidity_report_ja.tex').write_text(tex)
    for i in (1,2):
        with (OUT/('build-%d.log'%i)).open('w') as log:
            subprocess.run(['xelatex','-interaction=nonstopmode','-halt-on-error','cf_ep_rapidity_report_ja.tex'],cwd=str(OUT),stdout=log,stderr=subprocess.STDOUT,check=True)
    log=(OUT/'cf_ep_rapidity_report_ja.log').read_text()
    actual=re.search(r'Output written on .*? \((\d+) pages\)',log)
    if not actual:raise RuntimeError('No PDF page count')
    print('Requested sections:',pages,'PDF pages:',actual.group(1),'plot directory:',plotdir)
    for key in ['Missing character','Overfull','Underfull','Warning']:
        print(key,log.count(key))


if __name__=='__main__':main()
