# Λ質量fit：Poisson尤度Lによる再評価
現在の指定総合PDFは**L結果だけの20頁版**。[表示整理と旧結果の除外](poisson_L_only_report_20261102.md)を参照。以下の24頁版は直前の比較入り版の履歴であり、採用済みL数値自体は変更していない。

状態：完了。L＋非負定数背景fallbackで再評価し、総合PDF24頁を反映済み。

ユーザー追加指定：背景が負になる場合は、傾きを0に固定した非負の定数背景b≥0でL fitをやり直す。現行pol2では一次項・二次項をともに0へ固定する。主出力は`cf_poisson_likelihood_20261102_004/pdf/`。001–003は数値試行・定数fallback前の参照として保持する。

## 方針

- ユーザー指定の`L`を追加し、原カウントに`SLQNR0`を使用する。
- Gaussian＋pol2、mass fit範囲、初期値・パラメータ制限、最低観測数・共分散の判定を保持。
- fit範囲の空binを含むPoisson devianceと符号付きdeviance residualを計算する。
- Poissonの期待カウントが定義域外になる解は採用せず、その追加判定を区別して記録する。
- 連続fit範囲の背景最小値が負ならGaussian＋bへ再fitし、b≥0、自由parameter数4としてdeviance/NDFを評価する。fallbackのSE/ME識別を保存する。
- 5統合中央度700fitと9fine中央度1,260fitを同じ2021 P24iy 04BF ROOTから再抽出する。
- 前回の代表9点を固定してχ²版/L版を比較する。都合の良い点へ選び直さない。
- Method3/5/6を再計算し、共通規格化0.40–0.50で直接比較。0.30–0.40の復活も別途確認。
- 古い空bin除外pullを流用せず、空binも表示する。
- 総合PDFをL版へ更新し、χ²版PDFと図は退避・保存する。

## 作業

- [x] ROOT5で0×9bin＋1×1bin toyを確認（Lは0.1count/bin）。
- [x] L fit抽出と元カウントの一致・定義域・fit品質を検証。
- [x] CF、fine centrality純度、固定9点、全groupfit atlasを再作成。
- [x] χ²/Lの背景・純度・誤差・境界・devianceと空bin数を比較。
- [x] 規格化・共通ME純度誤差を再評価し図付きPDFへ反映。
- [x] Todo・README・日報・private snapshotを更新。

主出力：`cf_poisson_likelihood_20261102_004/pdf/`。
PID/Maker/mixing変更、production再投入、hadd、Git commit/pushは行わない。


---

## 完了結果（2026-11-02）

[直前の比較入り総合レポート（24頁）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_L_only_20261102_001/summary_4He_KFParticle_3p85.pdf)をL版へ更新し、[直前のχ²版（34頁）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_poisson_L_constant_20261102_004/summary_4He_KFParticle_3p85.pdf)を保存した。

### 最終的な背景の扱い

原カウントに`SLQNR0`でGaussian＋pol2をfitする。背景の最小値はfit範囲の端点と、範囲内にある二次式の頂点で評価する。背景が負なら、一次・二次係数をともに0に固定して**Gaussian＋非負定数b**でL fitをやり直す。bは0以上の範囲で推定し、Gaussianの3parameterは元の制限の範囲で自由とする。定数fitへ切り替えた場合は自由parameter数4、通常は6としてNDFを計算する。失敗した定数fitを負の元背景へ戻して採用することはしない。

多項式は中心化座標で同じ関数族を表し、係数・全共分散を元の表記へ変換して保存する。背景積分と共分散は安定な中心化基底で直接計算。初期背景は空binも含めたsideband平均（最低0.1count/binの初期値）、計算回数上限50,000。fit後のbの下限は0であり、0.1ではない。bの数値的上限max(1,Nfit)は全採用fitで非活性であることを確認した。

### 数値と代表例

1,960のmass-spectrum点を再評価。100点は最低観測数等によりfit未実行、1,860点のdevianceを独立照合。背景積分・共分散1,858点を検証。632点を定数背景で再fit（統合194、fine438）した。統合5中央度では700点中697点、fine9中央度では1,260点中1,158点でfull-window純度が有効。

| Method・中央度・k* [GeV/c] | χ²版P(SE) | L＋fallback P(SE) | 同じ点のP(ME) | SE/ME背景 |
|---|---:|---:|---:|---|
| M3 50–60%, 0.26–0.28 | 0.0104±0.4230 | 0.8531±0.0605 | 0.7605±0.0139 | 定数 / 二次 |
| M3 40–50%, 0.14–0.16 | 0.1850±0.3272 | 0.8813±0.0477 | 0.7849±0.0106 | 定数 / 定数 |
| M5 60–70%, 0.25–0.30 | 0.2510±0.2469 | 0.8180±0.0566 | 0.7107±0.00955 | 定数 / 二次 |
| M5 50–60%, 0.05–0.10 | 0.4571±0.1588 | 0.8297±0.0469 | 0.7212±0.0110 | 定数 / 定数 |

前回と同じ9例を固定して比較した。空bin除外に伴う大きな背景推定の偏りを避け、負の背景を非負定数で再評価できた。一方、SE/MEの一致は条件として課しておらず、差が残る点や高統計でのモデル残差もある。誤差の縮小には空bin情報に加え背景parameterの削減も寄与する。データが増えた、系統誤差まで解決した、という意味ではない。

規格化0.40–0.50、0.30–0.40の両方を再計算。比較の主区間は同じ0.40–0.50とし、以前不成立だった20–80%の旧区間も回復した。Method3/5/6の2分割・3分割CFを更新。747有効CF点と共分散、共有Method5 ME純度/errorを検証した。

### 図と検証

- [Method3 全fit atlas（125頁、0–0.5 GeV/c）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/atlas/atlas_method3.pdf)
- [Method5 全fit atlas（50頁、0–0.5 GeV/c）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/atlas/atlas_method5.pdf)
- [9中央度Method3：χ²/L比較](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/comparison/purity_fine_method3_comparison.pdf)
- [9中央度Method5：χ²/L比較](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/comparison/purity_fine_method5_comparison.pdf)
- [CF：χ²/L比較（同じ規格化、2頁）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/comparison/cf_objective_comparison.pdf)
- [数値・図の検証](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/validation.json)／[反映・元版保存](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/publication.json)／[設定・再現手順](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/README.md)

ROOT5 toyの9zero＋1oneはχ²で1count/bin、Lで0.0999948count/binを確認。12tests PASS。空binを落とす旧pullは使わず、全fit binの符号付きPoisson deviance residualを表示。全24頁を画像化し、表・タイトルの修正後はLaTeXのはみ出し警告なし。

b=0境界の局所共分散誤差は非対称な信頼区間の代わりにはならない。fineで1点がb≈0。モデル選択・mass model・ME→SE移行系統、未保存のevent/mixed-pair相関、背景fit/count共分散等は誤差に含めない。原ROOT・PID/Maker/mixingは変更せず、production再投入・hadd・Git commit/pushも行っていない。

最終保存：[個人snapshot 225件・hash照合](poisson_private_snapshot_20261102.json)／[最終確認](poisson_final_checks_20261102.json)。ROOT・バイナリ・大きな数値出力等は既存policyに従って除外し、元ファイルは保持した。commit/pushは行っていない。
