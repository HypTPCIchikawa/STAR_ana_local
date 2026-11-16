# EP分割・81 run較正計画の日本語PDF

**日付:** 2026-11-14 JST  
**対象:** Au+Au 3p85 / production_3p85GeV_fixedTarget_2021 / P24iy / Λ–⁴He  
**依頼:** EP分割の概念、今回の実装、未較正81 runで必要なパラメータ準備・検証を、学部生向けに図付き日本語PDFへまとめる。

## 成果物

- [PDF（A4縦・14ページ）](../../share/figure/20261114_ep_mixing_guide/ep_mixing_guide_ja.pdf)
- [TeX原稿](../../share/figure/20261114_ep_mixing_guide/ep_mixing_guide_ja.tex)
- [再生成・描画確認スクリプト](../../share/figure/20261114_ep_mixing_guide/build_pdf.py)
- [出典・SHA256・検証結果のmanifest](../../share/figure/20261114_ep_mixing_guide/manifest.json)
- [全ページ一覧](../../share/figure/20261114_ep_mixing_guide/contact_sheet.png)

PDF SHA256: `78ddb1634d48c1aefbf7be380ef12d7c597303ce50561d26ae96f147a1b10deb`、241850 bytes。

## 内容と判断の区別

1. EP/反応面、一次と二次の周期、SE/MEの向きの条件、flowを含む説明用モデル。
2. 実装済みEast EPD AB row1–8、nMIP≥0.3、重み上限2、normalized Q、run×較正centrality別recenter、20項shift、元Float精度のEP8分類。
3. 同じEP有効イベントでEP1/EP8を作る対照、cent9×Vz×EPのpool、最大20event・両方向mixing、raw SE/MEの合算保存。
4. 既存run22159051先頭10000eventの実EP図を再掲。valid9374、全8区分占有、SE一致の実装確認。全統計・較正合格・CF改善の証明とは区別。
5. 388 runと前期81 runの期間・coverage、Picoのgain補正と解析EP較正の違い、run list追加だけでは係数が生成されない構造。
6. 将来の共有候補／前期期間共通／run別係数の判断、recenter36平均とAB一次shift40係数、生成・適用条件の整合、独立検証と追加QA。

**81 runの解析・共有試験・新較正は保留のまま。** まず既存388 runでEP分割の有用性を確認し、有用と確認できたら81 run較正へ進む方針をPDFにも明記した。

## 検証

- 実装・較正の別々の読み取り調査後、原稿の物理内容を独立レビュー。修正必須なし。モデルの前提、残留平均の0との整合、EP1の各pool深さを明確化。
- `python3 build_pdf.py` PASS。XeLaTeX 2pass、14ページ、Missing character/Overfull/Underfull/Warningすべて0。
- Ghostscriptで全14ページを95dpi描画、全体一覧と主要図・最終参照ページを目視。日本語本文の抽出PASS。
- PDF圧縮streamも展開してDroidSansFallbackの日本語フォント資源、FontFile22、ToUnicode22を確認。
- 実EP図は元PDFと同一SHAで資料フォルダへ複製。資料単独での再生成に備えた。

## 出典と作業範囲

一般論はPoskanzer–Voloshin, Phys. Rev. C58 (1998), Sec. II [arXiv](https://arxiv.org/abs/nucl-ex/9805001)、Adams et al., The STAR Event Plane Detector, NIM A968 (2020) [arXiv](https://arxiv.org/abs/1912.05243)を照合。今回の値・仕様はsource/config/README/既存解析ノートから確認し、13出典ファイルのSHA256をmanifestへ保存。

解析source/config変更、新係数生成、Pico再読み取り、SUMS投入、追加メール送信は行っていない。進行中のhadd/通知のファイル・tmuxは変更していない。PDFに書いたhadd状態は12:40 JSTの観測時点を明記し、未出力groupを含む部分統計であることを保持した。

**mainconf / jobid / configlog:** 本PDF作成のための新規解析実行は無し。既存本番の記録は[EP1/EP8 hadd記録](pair_rapidity_hadd_20261114.md)を参照。
