# EP mixing区分別k* histogramの実装

## Session update — 20:29 JST / 2026-11-15：設定準備、C++読取承認待ち

ユーザーからEP8各planeのk*分布を保存できるよう解析コード修正の明示依頼。対象は2021 P24iy Λ–He4。新規本番jobは投入しない。

合意した出力は`hKstar{SE,ME}VsEpBinVsCent_lambda_he4{,_signal,_leftSB,_rightSB}`の8TH3D、軸はk*×実mixing EPclass×cent9。EP8は0–7、EP1対照は混合条件としてclass0だけを保存する。EP1で異なる物理Ψ1のeventを混ぜたMEを単一の物理planeに属すると解釈しない。既存SE/MEと同じ選択・weightを使い、各EP区分を足した際のclosureを検証する。

準備済み設定：

- mainconf: `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_kstarSplit_pairY10_covered25.yaml`
- histogram: `config/hist/hist_auau3p85_femtoLambda_he4_he3_control_epd8KstarPairY10_20261115_kf.yaml`
- analysis_info: `config/analysis/analysis_info_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_kstarSplit_pairY10_covered25.yaml`
- EP1検証fixture: `rootfile/femto_lambda_ep_kstar_validation_20261115_001/ep1_control_main.yaml`、`ep1_control_hist.yaml`
- 設定作成script/元設定SHA/新設定SHA: `rootfile/femto_lambda_ep_kstar_validation_20261115_001/prepare_configs.py`、`config_manifest.json`

元のhist定義・dataset・既存本番mainconfは変更せず、新設定を分離。較正済388runの入力catalogを継承、未較正81runは保留。histの存在によるopt-inと、全8本のTH3D型・既存対応histとのXZ軸一致・mixing設定とのY軸一致を初期化で検証する方針。FemtoConfigの新fieldは不要な構成を予定。

現時点の状態：C++本体は未変更。必要な限定ソース読取がroot/実装agent双方で自動承認レビューの再試行後もtimeoutし、同操作を停止。手動承認への一時切替をuserへ依頼済み。テスト担当は既に読めた情報から新規回帰の準備を続ける。build/実データ検証は未実施で、実装済みとは報告しない。

承認復旧後は、C++受け渡し・fill・init検証・EP境界metadataを実装し、batch整合build、EP8/EP1・全区分closure・ME両方向・invalid EP無充填の回帰と、可能なら小規模Pico検証を行う。入力ROOT/過去job snapshotを変更せず、メール・git pushは行わない。新mainconfに対するjobid/configlogは未生成。

追記：`tests/test_femto_epd_kstar_splits.C`と`tests/run_femto_epd_kstar_splits_test.C`をdraftとして保存。EP合算closure、両方向ME、EP1 class0、invalid class無部分充填、TH3のWrite後の保存を検証する設計。本体未実装のため未コンパイル・未実行であり、PASSではない。担当agentも追加作業を停止済み。

---

## Session update — 21:01 JST / 2026-11-15：ユーザー指示で中断

ユーザーが解析コードupdateの中断を明示。専用設定・テストdraftをそのまま保持し、C++本体の修正・build・追加検証・本番jobは実施しない。承認切替待ちによる停止ではなく、ユーザー方針による中断を最新状態とする。承認モードが変わっても自動再開しない。別作業として、保存済みROOTからcentrality別のpair rapidity PDFを作成する。
