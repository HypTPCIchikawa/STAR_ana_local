# Centrality別 Λ–He4 SE・ME・CF：ROOT直接PDF

[全体PDF：121ページ](centrality_pair_rapidity10_SE_ME_CF_root.pdf)。各centralityについて、inclusive、pair rapidityの区分1–10の順に11ページを配置しています。PDF bookmarkでも各条件へ移動できます。

各ページはROOTのCanvas Divide(2,3)です。左列EP1・右列EP8、上段SE＋scaled ME、中段CF全範囲、下段高k*拡大。タイトルにcentrality・cent9範囲・pair rapidity条件を記載しています。SE/ME総ペア数、規格化窓内のS_H/M_H、ME倍率、有効CF点数、空bin数も各列に表示します。これらは通常bin内のペア積分で、event数ではありません。

| Centrality | cent9 | 全体PDFのページ | 分冊PDF（各11ページ） |
|---|---|---|---|
| 0–20%（比較用合算） | 6–8 | 1–11 | [0–20%](by_centrality/lambda_he4_c0_20_SE_ME_CF_root.pdf) |
| 0–10%（比較用合算） | 7–8 | 12–22 | [0–10%](by_centrality/lambda_he4_c0_10_SE_ME_CF_root.pdf) |
| 0–5% | 8 | 23–33 | [0–5%](by_centrality/lambda_he4_c0_5_SE_ME_CF_root.pdf) |
| 5–10% | 7 | 34–44 | [5–10%](by_centrality/lambda_he4_c5_10_SE_ME_CF_root.pdf) |
| 10–20% | 6 | 45–55 | [10–20%](by_centrality/lambda_he4_c10_20_SE_ME_CF_root.pdf) |
| 20–30% | 5 | 56–66 | [20–30%](by_centrality/lambda_he4_c20_30_SE_ME_CF_root.pdf) |
| 30–40% | 4 | 67–77 | [30–40%](by_centrality/lambda_he4_c30_40_SE_ME_CF_root.pdf) |
| 40–50% | 3 | 78–88 | [40–50%](by_centrality/lambda_he4_c40_50_SE_ME_CF_root.pdf) |
| 50–60% | 2 | 89–99 | [50–60%](by_centrality/lambda_he4_c50_60_SE_ME_CF_root.pdf) |
| 60–70% | 1 | 100–110 | [60–70%](by_centrality/lambda_he4_c60_70_SE_ME_CF_root.pdf) |
| 70–80% | 0 | 111–121 | [70–80%](by_centrality/lambda_he4_c70_80_SE_ME_CF_root.pdf) |

pair YLabは[-1.5,0.5)の10等分、幅0.2です。各centralityの先頭ページはpair rapidityを制限しないinclusiveで、続く10ページが低いYから順に並びます。0–10%・0–20%は単独binと重複する比較用合算です。

規格化域は全条件で0.4 ≤ k* < 0.5 GeV/c、元200binを4binずつ加算した50binです。各EP・centrality・Y選択で別々のS_H/M_Hを用い、MEをSEの高さへ合わせています。合算centralityではSE・MEをそれぞれ加算してから規格化します。CFはsignal窓のraw CFで、背景差引・purity補正前です。表示範囲は0–1 GeV/c、高k*拡大は0.5–1 GeV/cです。

CFの誤差は共有規格化の誤差を含みます。event/pair再利用やEP1/EP8間の共分散は未評価です。SEまたはMEが空の点は省略し、規格化不能の条件には理由を表示します。とくに周辺centralityや一部Y区分の大きな誤差・空白を、CF=1の証拠と解釈しないでください。

2021 P24iyの既存hadd確定ROOT（両系列共通4023group）を使用しています。未較正81runは含みません。今回、production解析コードの更新・Pico再読・新規本番jobは行っていません。中断したEP別k*保存の実装と、この作図用postprocessingは別作業です。

検証は396投影closure、旧264ヒスト・132CF系列との完全一致、全484系列の独立数値検算、実際に描画した242列のペア数・有効点数の照合、PDF全121ページ描画です。[設定・来歴・SHA・検証manifest](canvas_manifest.json)、[PDF index](pdf_index.json)、[全centrality一覧画像](previews/centrality_overview.png)を保存しています。
