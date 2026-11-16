# 中央度・pair rapidityの結果を反映した日本語まとめ改訂

## Session update — 21:44 JST / 2026-11-15：33ページ版へ更新完了

ユーザー指定の[日本語まとめPDF](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/cf_ep_rapidity_report_ja.pdf)を22ページから33ページへ更新した。[旧22ページPDF](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/archive/20261115_initial22p/cf_ep_rapidity_report_ja.pdf)・旧TeX・manifest・組版sourceはSHA確認付きで保存。最新[詳細ROOT PDF121ページ](../../../share/figure/20261115_lambda_centrality_pair_canvases/centrality_pair_rapidity10_SE_ME_CF_root.pdf)と[11分冊一覧](../../../share/figure/20261115_lambda_centrality_pair_canvases/README.md)は保持。

## 追加した内容と結論

- 単独cent9全9区分と参照0–10/0–20%を説明し、SE総数・規格化窓・高k*窓の支持統計を明示。合算CFはSE/MEをそれぞれ足してから規格化し、CFの単純平均ではないことを記載。
- inclusiveの中央度別Ctail/T/GLS傾き、0–5%・20–30%のSE/ME/CF実例を追加。EP8 Ctailは0–5%0.89935±0.00421、10–20%0.92443±0.00681、20–30%0.94285±0.01027。0–30%はT<1と下降が残り、広い中央度の合算だけでは全体の残留を説明しない。
- EP1/EP8それぞれのpairY10×単独9中央度について、Ctail/T mapを追加。値±誤差を表示、NAは灰色、色域外と大きな誤差は表示上の注記で識別。元データの値・誤差・maskを変更していない。
- 固定pairY[-1.1,-0.7)の中央度依存を追加。中心側では傾きが小さくなるが水準に約2–3%残留、10–30%の1付近の結果は精度も含め限定標本として解釈。30–40%以降や端のpairYでの回復・再低下を断言しない。
- He4 Below/Aboveの中央度依存を追加。Belowでは10–20%に加え20–30%も現在の近似統計内で1・平坦に整合。Aboveは0–30%に残留。BelowのSE構成比も中央度で変わるため、inclusive変化を単一原因へ帰属しない。
- 規格化不能とTの可否の違い、周辺の少数count・Gaussian近似・event再利用相関未評価、固定normが物理的plateauとは未確認、k*上端1 GeV/cという限界を明記。結論と再現情報を新版に更新。

## 図・検証・来歴

追加ROOT vector図10本を`share/figure/20261115_lambda_he4_ep_rapidity_cf/centrality_update_plots/`へ保存。元180図は保持。図の504入力行の値/誤差/NAを最新JSONと独立照合し一致、ROOT5 exit0、全10図の描画/bbox検証を確認した。mapとsummaryに新しい物理選択やfitは加えていない。

独立科学レビューで新規中央度表・pairY傾向・He4二分・統計上の限界を確認。元studyの396closure、旧264hist/132CF完全一致、484系列/24,200点・共分散/542GLSの独立検算を継承。QA表には0–20%・EP8を明記し、過去の保留設定draftと今回のPDF後処理を区別した。

XeLaTeXの最終33ページはMissing character/Overfull/Underfull/Warningすべて0。全33ページを描画し一覧、mapとHe4形状の代表ページを目視。日本語text抽出・font埋込・ページ数・原ROOTのsize/mtimeNs/inode/device不変を確認後、指定PDFをatomic置換した。最終SHA256は`91b240a8221b756e20c8a9c6e4b954af376824e8c11f055d85ae9e4e5683b4fc`。

作業場所は`rootfile/femto_lambda_report_update_20261115_001/`。`centrality_scientific_summary.json/.md`、`summary_display.yaml`、`summary_plots.log`、`plot_input_provenance.json`、`summary_plot_manifest.json`、`centrality_pages.py`、`build_updated_report.py`、`verify_and_publish.py`、`publication_result.json`を保存。最終[report manifest](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/report_manifest.json)と[全33ページ一覧](../../../share/figure/20261115_lambda_he4_ep_rapidity_cf/contact_sheet.png)も更新。

## 中断方針と本番識別子

解析コードupdateはユーザー指示の中断を維持。既存EP実装準備の専用設定・未実行テストdraftは保持し、production C++は変更していない。本PDF更新のための新規本番mainconf/jobid/configlog、Pico再読・本番job・メール・pushは無し。未較正81runも保留。

入力は従来と同じ2021 P24iy、共通成功4023group・72424 Pico参照。mainconfは`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_covered25.yaml`、SUMS IDはEP1 `FD44A5EB1B0AA22C546F169E86482587` / EP8 `DA0DEB4347BD834F1096F8D8798F2133`、configlogは`job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_<jobid>.txt`。入力指紋は最新centrality studyの`input_manifest.json`を参照。
