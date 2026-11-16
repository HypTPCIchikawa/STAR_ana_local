# Λ–4He：SE/ME純度fit点検とMethod3・5・6のCF比較

2026-10-31（日本時間）。数値計算・全fit確認図・20頁report・成果物検証を完了した。


## 完成した図とレポート

- [Method3・5・6 比較とSE/ME fit点検レポート（20頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/method6_cf_purity_report_20261031.pdf)
- [CF重ね描き（中央度2分割・3分割、2頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/cf_methods_3_5_6.pdf)
- [SE/ME純度比較（規格化域まで表示、2頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/purity_se_me.pdf)
- [Method3 全fit atlas（100頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/atlas_method3.pdf)
- [Method5 全fit atlas（40頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/atlas_method5.pdf)
- [出力README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/README.md)／[全145PDFの頁数・hash](method6_outputs_20261031.json)／[最終検証](method6_validation_20261031.json)

## 対象とMethod6の定義

true=SE、mix=ME。2021 P24iy、SUMS `04BF2DAE47164B6F040E3A460E539289`、採用26,955/26,994 subjobs。10/30と同じ更新ROOTを読み取り専用で使用し、4He PID/TOF・Maker・mixingは変更しない。

ユーザー回答により、Method6は **Method5のME純度（50 MeV/c格子・滑動fit）とその誤差を、SE/MEの両方に共通適用**する。元5 MeV/c countへPを掛けてから20 MeV/cへ集約し、規格化を0.3–0.4 GeV/cで行う。

同じPと誤差を二つの独立推定として扱わず、同一の確率変数として共有する。純度cell内の規格化前SE/ME比ではPとその誤差が相殺する。50/20 MeV境界を跨ぐbinや、k*依存の重み付き規格化ではP依存が残る。ME→SE移行仮定や質量モデルの系統誤差は、コピーしたME fit統計誤差には含まれない。

## 実行条件と再現性

- 2分割：20–80%（cent9 0–5）、0–20%（6–8）。
- 3分割：0–10%（7–8）、10–20%（6）、20–60%（2–5）。
- CF 20 MeV/c、表示0–0.3 GeV/c、規格化0.3–0.4 GeV/c。
- Method3：各20 MeV/c binのSE/MEを独立fit。
- Method5：各50 MeV/c純度cellの前後2cellを含む最大250 MeV/c窓をSE/ME別にfit。
- atlasでは0≤k*<0.4の全fitを表示し、規格化域も含める。M3 100頁、M5 40頁、各頁はSE/MEのmass・model・BG・pullと正規化形状の6パネル。

1回のfit結果をJSONに保存し、CFと図の両方へ使った。描画時の再fitは行わない。10/30 Method3/5の500bin（欠測58binを含む）と一致し、最大値差・誤差差はともに約2.22×10⁻¹⁶。12 toy testsに加え、実データ198binで共有Pの局所相殺、642有効CF点と共分散を検証した。

[設定](settings/method6_settings.yaml)／[数値概要](method6_results_summary_20261031.json)／[出力ディレクトリ](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf)

## SEのestimationは正しいか

**Method3の最初のk* binはSEの低統計が明確な問題**。例えば20–80%の0–0.02 GeV/cはmass全域でもSE観測3件でfitできない。0–20%でも22件しかなく、数値収束しても背景が非物理的になり純度を採用できない。

**Method5のSE/ME差を、SE低統計だけでは説明しきれない**。中心cellが少数でも、滑動fitは広い窓の統計を使う。0–20%の最初の50 MeV/c cellでは、SEのfit窓0–0.15に11,633件があり、P_SE=0.916347±0.001565、P_ME=0.920864±0.000249。χ²/NDFはSE約2.66、ME約43.88で、MEにもモデル不適合がある。

さらに同じ群の0.30–0.35 cellでは、SE fit窓270,578件、ME約1,001万件に対し、χ²/NDFはSE約38.54、ME約1329.24。統計誤差が小さいことはモデルや純度の正しさを保証しない。SEだけを不良、MEを正解と認定しない。

下表は表示・規格化0–0.4内の、数値的に利用できる純度件数とfit適合度。利用可はモデル妥当性の認定ではない。

| 中央度 | 方法 | SE純度可 | ME純度可 | median χ²/NDF SE | ME |
|---|---|---:|---:|---:|---:|
| 20-80% | method3 | 18/20 | 20/20 | 1.092 | 4.993 |
| 20-80% | method5 | 8/8 | 7/8 | 3.536 | 60.174 |
| 0-20% | method3 | 19/20 | 20/20 | 2.333 | 33.133 |
| 0-20% | method5 | 8/8 | 8/8 | 16.701 | 504.469 |
| 0-10% | method3 | 19/20 | 20/20 | 2.261 | 27.734 |
| 0-10% | method5 | 8/8 | 8/8 | 13.449 | 421.239 |
| 10-20% | method3 | 18/20 | 20/20 | 1.263 | 6.832 |
| 10-20% | method5 | 8/8 | 8/8 | 4.593 | 88.416 |
| 20-60% | method3 | 18/20 | 20/20 | 1.092 | 5.012 |
| 20-60% | method5 | 8/8 | 8/8 | 4.085 | 60.146 |

## SE/ME純度一致の仮定について

SE内で信号と背景の相関が異なれば、十分な統計でもSE/MEの純度は異なり得る。ME純度をP、信号／背景の相関をC_s/C_bと書けば、SE純度は `C_s P / [C_s P + C_b (1−P)]` となる。差だけでfit不良とは断定できないため、質量形状・model・統計の点検とMethod6比較を分けて扱う。

Method6はユーザー指定のME→SE移行仮定の比較であり、独立なSE背景推定や最終採用の認定ではない。誤差が小さくなることも仮定の正しさを保証しない。異なる方式間の差の有意性は計算していない。

## CF比較と欠測

| 中央度 | Method3積分CF | Method5積分CF | Method6積分CF |
|---|---:|---:|---:|
| 20-80% | NA | NA | NA |
| 0-20% | NA | 1.024121 ± 0.008792 | 1.026155 ± 0.008551 |
| 0-10% | NA | 1.024593 ± 0.009838 | 1.027566 ± 0.009589 |
| 10-20% | NA | 1.044588 ± 0.019951 | 1.030500 ± 0.019098 |
| 20-60% | NA | 1.101549 ± 0.026633 | 1.077379 ± 0.024990 |

積分は0≤k*<0.2 GeV/cのcount積分から求め、0.3–0.4で規格化する。点の単純平均ではない。Method3は低k*欠測により全群で積分がNAだが、有効なCF点は図に残る。

20–80%はMEの0.30–0.35 purity fitが不成立のため、Method5・6とも規格化不可でCF全体をNAとする。MEを共有するだけでは回復しない。norm域を縮めたり失敗fitを0へ置き換えたりしない。

## 検証・保存

入力fingerprint、manifest/runmeta/runtime bundle、config/source/copy hashを照合。元ROOT・旧reportと数値は保持する。未保存のfit-background/count・異なる滑動fit間の背景共分散・event/mixed-pair再利用相関は完全には評価できず、誤差は一次近似。

図の余白は元のvector PDFを保存した上でCropBoxで調整し、data・bin・fitは変更しない。評価コードと小さなPDF/図は個人サブモジュールへ保存し、共有Gitへpushしない。

完了時刻：2026-10-31T13:02:48.744193+09:00。

個人保存：[今回のsnapshot一覧](method6_private_snapshot_20261031.json)。157件を保存し、ROOT/build・大きな数値出力・1 MiBを超える全fit集PDFは方針に従い除外した。全fit集の各1頁PDFと20頁reportは小さいため保存した。共有Gitへのcommit/pushは行っていない。
