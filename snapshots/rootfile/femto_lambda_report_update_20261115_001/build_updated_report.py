#!/usr/bin/env python3
"""Build a Japanese report that includes the centrality-resolved CF study."""
from pathlib import Path
import importlib.util
import json
import re
import subprocess
import sys
import yaml

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
STAGE=HERE/'pdf_build'
sys.path.insert(0,str(HERE))
from centrality_pages import Pages, pm


def header(title,kind):
    return '\\phead{00}{'+title+'}{'+kind+'}\n'


def new_overview(pages):
    t=[header('EP・中央度・rapidityで高k* CFを検証する','中央度別結果を反映した改訂版')]
    t.append(r'\lead{EP8で改善するが、全体の残留は解決していない。}')
    t.append(r'同じ規格化$H=[0.4,0.5)$ GeV/$c$で、signal raw CFを比較しました。$C_{\rm tail}$は$0.8\leq k^*<1$ GeV/$c$のcountを積分してから作った比です。')
    t.append(r'\begin{center}\small\begin{tabular}{P{6.4cm}cc}\toprule 条件 & EP1 & EP8\\\midrule')
    for label,cent,sel in [
        (r'0--20\% inclusive','c0_20','inclusive'),
        (r'0--5\% inclusive','c0_5','inclusive'),
        (r'20--30\% inclusive','c20_30','inclusive'),
        (r'0--20\%：pairY $[-1.1,-0.7)$','c0_20','pair_merge02'),
        (r'20--30\%：He4 Below','c20_30','he4_below')]:
        t.append(label+' & '+pages.tail(cent,sel,'EP1')+' & '+pages.tail(cent,sel,'EP8')+r'\\')
    t.append(r'\bottomrule\end{tabular}\end{center}')
    t.append(r'\sub{中央度だけを狭めても、下降は残る}')
    t.append(r'0--20\%全体は0.866から0.906へ改善。それでも約9.4\%低く、単独0--5\%や10--20\%、20--30\%にも下降が残ります。広い中央度を合算したことだけでは説明できません。')
    t.append(r'\sub{ラピディティで残留と精度の両方が変わる}')
    t.append(r'固定pairY $[-1.1,-0.7)$では中心側の傾きが小さく、He4 Belowの10--20\%・20--30\%は限定標本として1・平坦と整合します。周辺中央度や端のYでは低統計が大きく、全体の回復には外挿しません。')
    t.append(r'\warmnote{評価は保存ヒストの上端までの\textbf{0.5--1 GeV/$c$}です。固定した規格化域が物理的plateauとは未確認であり、1 GeV/$c$超の漸近収束や残留原因の特定は行っていません。}')
    t.append(r'\figcap{誤差は共有規格化を含むhistogram統計近似。event/pair再利用、EP1/EP8間の共分散、系統誤差は未評価です。解析コードupdateは中断を維持しています。}')
    return '\n'.join(t)+'\n'


def new_conclusion():
    return header('結論：中央度を加えて何が分かったか','更新した判断')+r'''
\lead{EP8の改善、残るY依存、低統計の限界を分ける。}
\sub{1. 全体の問題は未解決}
0--20\% inclusiveの$C_{\rm tail}$は0.866から0.906へ上がりました。単独の0--5\%、5--10\%、10--20\%、20--30\%でも、EP8に1未満と下降が残ります。中央度の幅を狭めるだけで高k*の形状差が解消する結果ではありません。
\sub{2. 同じ中央度の中にもpairY依存がある}
0--30\%では、pairYが$[-0.7,-0.5)$から$[-0.5,-0.3)$、$[-0.3,-0.1)$へ大きくなると水準が下がる傾向が見えます。固定$[-1.1,-0.7)$は傾きが小さくなりますが、中心側には約2--3\%の水準残留があります。
\sub{3. He4 Belowの整合性は限定標本の結果}
10--20\%に加えて20--30\%も、He4 Belowでは現在の近似統計内で1・平坦と整合します。一方、Aboveは0--30\%に残留が続きます。中央度でY標本の構成も変わるため、inclusiveの変化を単一の原因へ帰属しません。
\sub{4. 周辺中央度の色や誤差棒だけで回復を決めない}
規格化窓や高k*窓のcountが減り、NAや大きな誤差が増えます。50--60\%以降の一部の低い中心値を再低下と断言せず、1をまたぐ誤差だけを回復の証拠にもしません。
\boxnote{今回のY分割はSE/ME共通のペア選択です。Y別に独立したME poolを作る効果や、EP8各plane別CFは検証していません。EP別k*保存の解析コードupdateはユーザー指示で中断中です。新しいPico処理・本番job・未較正81 runの較正を行っていません。}
\sub{今後検討する項目}
固定norm域の依存性、残留の大きい／小さい条件の運動学・受容・ME形状の比較、高k*範囲の拡張、event再利用相関の評価が候補です。これらは未実施であり、現在の図からMEだけを原因と断定しません。
'''


def new_appendix():
    return header('最新版の再現と詳細図の参照先','来歴・検証・中断事項')+r'''
\lead{入力ROOTと旧版を保持し、後処理だけを更新。}
2021 P24iy、較正対象388 runの共通成功4,023 group、72,424 Pico参照を継承。未出力35 groupと未較正81 runは含みません。本番識別子はEP1 \texttt{FD44A5EB1B0AA22C546F169E86482587}、EP8 \texttt{DA0DEB4347BD834F1096F8D8798F2133}です。
\sub{元の0--20\%中心の結果}
{\small\texttt{rootfile/femto\_lambda\_ep\_rapidity\_cf\_20261115\_001/}}
132 CF系列、264 raw投影、180 vector図を保持しています。旧22ページPDF・TeX・manifestは本PDFのdirectory内の\texttt{archive/20261115\_initial22p/}に保存しました。
\sub{中央度を拡張した最新study}
{\small\texttt{rootfile/femto\_lambda\_cent\_pair\_canvas\_20261115\_001/}}
\texttt{study.yaml}に11slice・norm・rebin、\texttt{analysis\_results.json}に484系列、\texttt{metrics.csv}と\texttt{cf\_covariances.npz}に定量と共分散を保存。\texttt{independent\_validation.json}に独立再計算、\texttt{cf\_histograms.root}にSE/ME/scaled ME/CF graphを保存しています。
\sub{詳細ROOT PDF}
{\small\texttt{share/figure/20261115\_lambda\_centrality\_pair\_canvases/}}
\texttt{centrality\_pair\_rapidity10\_SE\_ME\_CF\_root.pdf}は121ページ。各centralityの11ページ分冊・README索引・\texttt{pdf\_index.json}も同directoryにあります。
\sub{今回のまとめ更新}
{\small\texttt{rootfile/femto\_lambda\_report\_update\_20261115\_001/}}
中央度の科学レビュー、図の表示YAML・入力SHA・作図log、新組版sourceを保存。追加vector図は本PDFのdirectory内の\texttt{centrality\_update\_plots/}、最終来歴は\texttt{report\_manifest.json}を参照してください。
\sub{検証と限界}
396投影closure、旧264ヒスト・132CF系列の完全一致、全484系列・24,200点と共分散・542 GLSの独立検算を継承。元merge ROOTはsize/mtimeNs/inode/deviceで不変を確認しています。元ROOT全file SHA、event再利用共分散、0--80\%の独立投影比較は未実施です。
\figcap{ROOT5.34/38・SL24y wrapper、system Python3 isolated mode/NumPy1.23.5、XeLaTeXを使用。本PDF更新のための本番mainconf/jobid/configlog・メールは新規作成していません。EP実装準備の専用設定・未実行テストdraftは、解析コードupdateの中断を維持して保持しています。}
'''


def build():
    settings=yaml.safe_load((HERE/'report_settings.yaml').read_text())
    settings['plot_files']={
        'inclusive_tail':'inclusive_ctail_vs_centrality.pdf','inclusive_shape':'inclusive_shape_vs_centrality.pdf',
        'ep1_tail_map':'EP1_pairY_ctail_map.pdf','ep8_tail_map':'EP8_pairY_ctail_map.pdf',
        'ep1_shape_map':'EP1_pairY_shape_map.pdf','ep8_shape_map':'EP8_pairY_shape_map.pdf',
        'merge_tail':'fixedmerge02_ctail_vs_centrality.pdf','merge_shape':'fixedmerge02_shape_vs_centrality.pdf',
        'he4_tail':'he4_splits_ctail_vs_centrality.pdf','he4_shape':'he4_splits_shape_vs_centrality.pdf'}
    (HERE/'report_settings.yaml').write_text(yaml.safe_dump(settings,sort_keys=False,allow_unicode=True))
    source=ROOT/settings['original_study']/'build_report.py'
    spec=importlib.util.spec_from_file_location('original_report',source)
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    tex,original_pages,_=base.build_tex();assert original_pages==22
    sections=re.split(r'(?=\\phead\{\d+\})',tex);preamble=sections[0];old={}
    for section in sections[1:]:
        number=int(re.match(r'\\phead\{(\d+)\}',section).group(1));old[number]=section
    assert len(old)==22
    pages=Pages(ROOT,settings);additional=pages.build();assert len(additional)==13
    old[1]=new_overview(pages)
    old[2]=old[2].replace(
        r'主結果は0--20\%（cent9 6--8）。cross-checkは0--10\%（7--8）と10--20\%（6）。',
        r'0--20\%を参照に、単独cent9全9区分と0--10\%合算を比較。単独中央度の対応・支持統計は後の追加章に示します。')
    old[19]=old[19].replace(r'\sub{閉包は108件すべて合格}',r'\sub{中央度拡張後の396閉包も合格}')
    old[19]=old[19].replace('\n','\n'+r'\figcap{以下の表は0--20\%・EP8の支持統計です。全中央度のpairY別SE$_H$は、前の統計表に示しています。}'+'\n',1)
    old[19]=old[19].replace(r'今回、$0\leq k^*<1$のSE/MEでpairY範囲外のcountはありません。',
                          r'単独9中央度と合算2条件の$0\leq k^*<1$で、SE/MEのpairY範囲外countは0でした。')
    old[19]=old[19].replace(r'0--20\%＝0--10\%＋10--20\%も独立に確認しています。',
                          r'0--10\%＝0--5\%＋5--10\%、0--20\%＝0--5\%＋5--10\%＋10--20\%も内容・分散を独立に確認。旧264ヒスト・132CF系列を再現しました。')
    old[21]=new_conclusion();old[22]=new_appendix()
    ordered=[old[i] for i in range(1,15)]+additional+[old[i] for i in range(17,23)]
    assert len(ordered)==33
    for i,section in enumerate(ordered,1):
        ordered[i-1]=re.sub(r'\\phead\{\d+\}',lambda _:r'\phead{%02d}'%i,section,count=1)
    tex=preamble+'\n'.join(ordered)+r'\label{reportend}\end{document}'+'\n'
    tex=tex.replace(r'\vspace{2mm}\small #1',r'\vspace{2mm}\small\raggedright #1')
    tex=tex.replace('EP・rapidity比較 / 2026年11月15日 JST','EP・中央度・rapidity / 2026年11月15日 JST')
    old_out=ROOT/'share/figure/20261115_lambda_he4_ep_rapidity_cf'
    tex=re.sub(r'\{(plots_cropped/[^}]+)\}',lambda match:'{'+str(old_out/match.group(1))+'}',tex)
    tex=tex.replace(r'\label{reportend}\end{document}\label{reportend}\end{document}',r'\label{reportend}\end{document}')
    assert len(re.findall(r'\\phead\{\d+\}',tex))==33
    return tex,settings


def main():
    STAGE.mkdir(exist_ok=True);tex,settings=build()
    (STAGE/'cf_ep_rapidity_report_ja.tex').write_text(tex)
    for i in (1,2):
        with (STAGE/('build-%d.log'%i)).open('w') as log:
            subprocess.run(['xelatex','-interaction=nonstopmode','-halt-on-error','cf_ep_rapidity_report_ja.tex'],
                           cwd=str(STAGE),stdout=log,stderr=subprocess.STDOUT,check=True)
    log=(STAGE/'build-2.log').read_text();match=re.search(r'Output written on .*? \((\d+) pages\)',log)
    if not match or int(match.group(1))!=33:raise RuntimeError('Unexpected PDF page count')
    warnings={key:log.count(key) for key in ('Missing character','Overfull','Underfull','Warning')}
    if any(warnings.values()):raise RuntimeError('Typesetting warnings: '+str(warnings))
    (HERE/'build_summary.json').write_text(json.dumps(dict(pages=33,warnings=warnings,staged_pdf=str(STAGE/'cf_ep_rapidity_report_ja.pdf')),indent=2)+'\n')
    print('UPDATED_REPORT_BUILT pages=33 warnings=0 staged_only=YES')


if __name__=='__main__':main()
