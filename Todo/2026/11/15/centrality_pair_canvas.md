# centrality別pair rapidity ROOT直接PDF

- [x] 解析コードupdateはユーザー指示で中断し、既存設定・テストdraftを保持する。
- [x] cent9対応を確認：0–5,5–10,10–20,20–30,30–40,40–50,50–60,60–70,70–80%の9区分、比較用0–10/0–20を加えた11区分。
- [x] 保存済みEP1/EP8 ROOTから11centralityのSE/MEを投影し、各条件の[0.4,0.5)規格化・CFを計算（484系列、396closure PASS）。
- [x] pairY10区分とinclusiveについて、各centralityのSE+scaled ME/CF full/highをROOT Divide(2,3)で出力。全121ページ＋centrality別11分冊各11ページ。
- [x] 全121ページ描画・旧264hist/132CF完全一致・全CF独立検算・描画242列のcounts一致・原ROOT指紋不変を確認し、PDF/ROOT/manifestを記録。

21:01 JST。production source変更・Pico再読・本番job無し。未較正81runは保留。

21:12 JST完了。[PDF・分冊一覧・検証記録](../../../../analysisnote/20261115/centrality_pair_canvas.md)。解析コードupdateの中断を維持。
