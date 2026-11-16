"""Japanese report pages generated from the measured centrality CF results."""
from pathlib import Path
import json
import math


def pm(value, error):
    if value is None or error is None:
        return '不可'
    digits=max(0,min(4,1-int(math.floor(math.log10(error))))) if error>0 else 4
    return '$%.*f\\pm%.*f$'%(digits,value,digits,error)


class Pages:
    def __init__(self, root, settings):
        self.root=root;self.settings=settings
        self.data=json.loads((root/settings['centrality_study']/'analysis_results.json').read_text())
        self.singles=settings['single_centralities'];self.sections=[];self.chunks=[]
        self.plotnames=settings['plot_files']

    def r(self, cent, selection='inclusive', dataset='EP8'):
        return self.data['curves']['__'.join((dataset,cent,selection))]

    def tail(self, cent, selection='inclusive', dataset='EP8'):
        r=self.r(cent,selection,dataset)['windows']['k080_100'];return pm(r['value'],r['error'])

    def shape(self, cent, selection='inclusive', dataset='EP8'):
        r=self.r(cent,selection,dataset)['shape_ratio'];return pm(r['value'],r['error'])

    def put(self, text):self.chunks.append(text+'\n')

    def page(self, title, kind='中央度の追加検証'):
        if self.chunks:self.sections.append(''.join(self.chunks))
        self.chunks=['\\phead{00}{'+title+'}{'+kind+'}\n']

    def table(self, columns, heading, rows, size='small'):
        self.put('\\begin{center}\\'+size+'\\begin{tabular}{'+columns+'}\\toprule\n'+heading+r'\\\midrule')
        for row in rows:self.put(' & '.join(row)+r'\\')
        self.put(r'\bottomrule\end{tabular}\end{center}')

    def image(self, path, width='0.94\\linewidth', page=None):
        path=Path(path)
        if not path.is_file():raise FileNotFoundError(path)
        options='pagebox=cropbox,width='+width
        if page is not None:options='page='+str(page)+','+options
        self.put(r'\begin{center}{\normalfont\special{pdf:literal direct q 0 Tr}'
                 +'\\includegraphics['+options+']{'+str(path)+'}'
                 +r'\special{pdf:literal direct Q}}\end{center}')

    def plot(self, key, width='0.94\\linewidth'):
        self.image(self.root/self.settings['plots_directory']/self.plotnames[key],width)

    def finish(self):
        if self.chunks:self.sections.append(''.join(self.chunks));self.chunks=[]
        return self.sections

    def build(self):
        self.page('中央度を9区分に広げて確かめる','更新範囲・統計量')
        self.put(r'\lead{0--20\%の合算結果と、単独cent9の結果を比較する。}')
        self.put(r'保存済みROOTを0--5、5--10、以後10\%刻みで70--80\%まで分割しました。従来の0--10\%・0--20\%も参照として残します。単独9区分と合算2区分は重複するため、11標本を独立とは扱いません。')
        rows=[]
        for c in self.singles:
            r=self.r(c['id']);rows.append([c['label'].replace('%',r'\%'),str(c['min']),format(r['totalSE'],',.0f'),format(r['normSE'],',.0f'),format(r['windows']['k080_100']['SE'],',.0f')])
        self.table('lrrrr',r'中央度 & cent9 & EP8 SE総数 & SE$_H$ & SE$_{[0.8,1)}$',rows)
        self.put(r'\figcap{countは$0\leq k^*<1$の通常binに入ったペア積分。event数ではありません。周辺中央度では、特に規格化窓と高k*窓の統計が減少します。}')
        self.put(r'\sub{中央度を合算する手順}')
        self.put(r'各cent9内でMEを作った後、$S=\sum_c S_c$、$M=\sum_c M_c$を作り、その比を規格化します。cent9をまたぐevent mixingは追加しません。合算後の比は、各CFの単純平均とは重みが異なります。')
        self.put(r'\boxnote{今回増やしたのは保存済みヒストの投影条件です。解析コードupdateはユーザー指示で中断しており、新しいPico処理・本番job・81 runの較正は行っていません。}')

        self.page('全体CFの水準には中央度依存がある','inclusive / 単独9区分')
        self.plot('inclusive_tail')
        rows=[[c['label'].replace('%',r'\%'),self.tail(c['id'],dataset='EP1'),self.tail(c['id'])] for c in self.singles]
        self.table('lcc',r'中央度 & $C_{\rm tail}$ EP1 & $C_{\rm tail}$ EP8',rows)
        self.put(r'\figcap{横軸は中央度。大きい側ほど周辺衝突です。EP8は0--5\%で約0.899、10--20\%で約0.924、20--30\%で約0.943。30--40\%まで中心値は1へ近づきますが、1への普遍的な回復とは結論しません。}')
        self.put(r'\warmnote{50--60\%以降の精度低下を含め、全域で単調な中央度依存が確定したという意味ではありません。70--80\%は高k*窓のSE=0で評価不可です。}')

        self.page('水準が上がっても、下降は残っている','inclusive / 形状と傾き')
        self.plot('inclusive_shape')
        rows=[]
        for c in self.singles:
            f=self.r(c['id'])['fits']['high']
            rows.append([c['label'].replace('%',r'\%'),self.shape(c['id']),pm(f.get('slope'),f.get('slope_error'))])
        self.table('lcc',r'中央度 & EP8 $T$ & EP8傾き $(\mathrm{GeV}/c)^{-1}$',rows)
        self.put(r'\figcap{$T$は高い側・低い側の2窓の積分$S/M$の比です。0--30\%では$T<1$と負の傾きが残ります。単独0--5\%でも低下するため、広い0--20\%を合算したことだけでは残留を説明できません。}')
        self.put(r'\boxnote{周辺中央度の大きな誤差・少数点から、平坦化や再低下を確定しません。傾きは共分散を含む0.5--1 GeV/$c$の診断量です。}')

        index=json.loads((self.root/'share/figure/20261115_lambda_centrality_pair_canvases/pdf_index.json').read_text())
        main_canvas=self.root/index['combined_pdf']
        lookup={x['centrality']['id']:x for x in index['centralities']}
        for cent,title,caption in [
            ('c0_5','SEとMEの形状差：最中央0--5\%',r'0--5\%ではSE/MEの差と高k*の下降が、単独cent9内でも残っています。'),
            ('c20_30','SEとMEの形状差：20--30\%',r'20--30\%ではEP8で差が小さくなりますが、全体CFの高k*下降が残ります。')]:
            self.page(title,'ROOT出力の代表例 / inclusive')
            self.image(main_canvas,'0.98\\linewidth',lookup[cent]['combined_first_page'])
            self.put(r'\figcap{'+caption+r'左EP1・右EP8、上SEとscaled ME、中CF全体、下高k*拡大。各列のS/M総ペア数と規格化窓のcountも併記しています。}')
            self.put(r'\boxnote{EP8の$C_{\rm tail}$='+self.tail(cent)+r'、$T$='+self.shape(cent)+r'。規格化域は共通の$[0.4,0.5)$ GeV/$c$です。}')

        for ds in ('EP8','EP1'):
            for kind,name in [('tail',r'$C_{\rm tail}$'),('shape',r'$T$')]:
                self.page('pairYと中央度の対応：'+ds+' / '+('CF水準' if kind=='tail' else '形状比'),'pairY10 × 単独9中央度')
                self.plot(ds.lower()+'_'+kind+'_map','0.98\\linewidth')
                self.put(r'\figcap{各cellは1つの中央度・pair rapidity条件で、値と統計近似誤差を併記しています。灰色NAは計算不可。色域外は端色と斜線で示し、表示範囲の全幅以上の誤差には!を付けています。これらは表示上の注記で、追加の選択cutではありません。}')
                if kind=='tail':
                    self.put(r'\sub{同じ中央度の中でも、Yが変わるとCFが変わる}')
                    self.put(r'0--30\%では、pairYの$[-0.7,-0.5)$、$[-0.5,-0.3)$、$[-0.3,-0.1)$と大きい側へ進むと、CF水準が低くなる傾向が見えます。中央度だけを分けても、このY依存は消えません。')
                else:
                    self.put(r'\sub{水準と形状を両方見て判断する}')
                    self.put(r'$T=1$は2窓の積分$S/M$が等しいことを表し、窓内まで平坦である証明ではありません。Yが大きい側には下降が残る一方、$[-1.1,-0.7)$付近では形状差が小さくなります。')
                self.put(r'\warmnote{40--50\%以降や端のY区分には低統計・NAが多く、色だけから中央度の連続的な変化や回復を読み取りません。EP1/EP8のmapは同じ色域を使用しています。}')

        self.page('固定pairY区間の「ほぼ平坦」はどこまで続くか','pairY [-1.1,-0.7) / 事前固定merge')
        self.plot('merge_tail','0.87\\linewidth')
        self.plot('merge_shape','0.87\\linewidth')
        self.put(r'\figcap{既存の隣接2bin mergeを全中央度で同じように適用しました。0--5\%は$C_{\rm tail}$='+self.tail('c0_5','pair_merge02')+'、$T$='+self.shape('c0_5','pair_merge02')+r'。中心側には約2--3\%の水準残留があっても、傾きは小さくなります。}')
        self.put(r'\figcap{10--20\%は'+self.tail('c10_20','pair_merge02')+'、20--30\%は'+self.tail('c20_30','pair_merge02')+r'。統計誤差も増えるため、0--20\%の改善を全中央度へ外挿しません。30--40\%以降は支持統計がさらに減少します。}')

        self.page('He4自身のY選択：中央度ごとの水準','He4 Below / Above')
        self.plot('he4_tail','0.95\\linewidth')
        rows=[[c['label'].replace('%',r'\%'),self.tail(c['id'],'he4_below'),self.tail(c['id'],'he4_above')] for c in self.singles[:5]]
        self.table('lcc',r'中央度 & EP8 Below & EP8 Above',rows)
        self.put(r'\figcap{Belowは$y_{\rm He,Lab}<-0.5$、Aboveは$>-0.5$。Belowの10--20\%・20--30\%は1付近ですが、Aboveでは0--30\%に残留が続きます。pairY選択との同時条件を課した結果ではありません。}')

        self.page('He4自身のY選択：形状と標本構成','限定標本の整合性')
        self.plot('he4_shape','0.85\\linewidth')
        rows=[[c['label'].replace('%',r'\%'),self.shape(c['id'],'he4_below'),self.shape(c['id'],'he4_above')] for c in self.singles[:4]]
        self.table('lcc',r'中央度 & EP8 Below $T$ & EP8 Above $T$',rows,'footnotesize')
        self.put(r'\figcap{20--30\% Belowでは$T$='+self.shape('c20_30','he4_below')+r'、傾き='+pm(self.r('c20_30','he4_below')['fits']['high']['slope'],self.r('c20_30','he4_below')['fits']['high']['slope_error'])+r'。10--20\%に加えて20--30\%も、現在の近似統計内で1・平坦と整合する限定条件です。}')
        fractions=[]
        for c in self.singles[:4]:
            den=self.r(c['id'])['totalSE'];fraction=self.r(c['id'],'he4_below')['totalSE']/den
            fractions.append(c['label'].replace('%',r'\%')+'：'+('%.1f'% (100*fraction))+r'\%')
        self.put(r'\boxnote{inclusive SEに占めるBelowの比率も中央度で変わります（'+ '、'.join(fractions)+r'）。inclusive CFの中央度依存を、EPだけ・Yだけなど単一の原因に帰属できません。}')

        self.page('どの区分に十分な規格化統計があるか',r'EP8 / pairY10の$S_H$')
        rows=[]
        for c in self.singles:
            rows.append([c['label'].replace('%',r'\%')]+[format(self.r(c['id'],'pair_bin%02d'%i)['normSE'],',.0f') for i in range(1,11)])
        self.table('l'+'r'*10,r'中央度 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10',rows,'footnotesize')
        self.put(r'\figcap{列番号はpairYを$[-1.5,0.5)$で10等分したものです。区分$i$は$[-1.5+0.2(i-1),-1.5+0.2i)$。SE$_H$は規格化窓$[0.4,0.5)$内のSE積分です。ME$_H$と高k*窓の統計も、CFを判断するために必要です。}')
        self.put(r'\sub{値が出せない理由を区別する}')
        self.put(r'50--60\%の固定pairY mergeではSE$_H=0$のためCFが計算できません。一方、規格化が相殺される$T$は大きな誤差付きで計算できる場合があります。CFと$T$の可否は同一ではありません。')
        self.put(r'\sub{周辺中央度では誤差近似にも注意する}')
        self.put(r'70--80\% inclusiveはSE$_H=2$、$[0.8,1)$のSE=0です。規格化窓に少数ペアがあるだけでは、高k*のCFや平坦性を判断できません。少数countのGaussian誤差やGLSは診断量として扱います。')
        self.put(r'\boxnote{中心値が1に近い、あるいは誤差棒が1をまたぐ、ということだけを一般的な回復の証拠にしません。支持count、形状比、実際の図を合わせて読みます。}')
        return self.finish()
