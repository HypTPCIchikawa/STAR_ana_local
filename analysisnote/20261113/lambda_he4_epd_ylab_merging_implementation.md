# Λ–⁴He：EP・YLab・track-merging解析コードの反映

記録：2026-11-13 JST（study名は作業開始UTC日の20261112）。対象は2021年P24iy、production_3p85GeV_fixedTarget_2021。新規commit・SUMS投入・mail送信なし。

## 実装結果

- Xialei方式のEast EPD AB（row1–8）一次Ψ1をrecenter＋20項shiftで補正し、元Float精度の8区分をcent9×Vzのmixing keyへ追加した。EP1対照も同じEP有効eventだけを使う。TPCの補助Ψ2とは別で、代用しない。イベントごとのtile samplingは通常解析の乱数streamを消費しない。元の乱数履歴がないため、元treeの角度とのbitwise一致とは主張しない。
- 各pairの物理⁴He自身のP4から求めたYLabで、strict `YLab < -0.5` / `YLab > -0.5` のraw k* SE/MEを追加した。全質量・signal・左右sidebandの各TH1とk*×cent9、計32ヒスト。ME両方向とも⁴He側の条件を使う。He3AsHe4は除外。等号・非有限値は追加diagnosticだけに入れず、inclusive物理pairは維持する。同eventに両側のHe4があれば両側へ寄与するので、排他的なevent分類ではない。
- `FillMerging`の旧primary pMom/pPt gate、誤った曲率符号・固定磁場・Heの二重Z補正・asin clampを除き、実event磁場[kG]とglobal rigidityに基づくSTAR helixで半径1.4mへの到達位置を計算する。非負の最短path・半径残差0.001cm以内・|z|≤210cmを要求し、到達しないtrackは除く。Δφ*は位置φ差、Δηはglobal運動量η差。protonに加えてpion図を新profileへ追加。bareと_EventFieldキーは修正後の同じ定義の互換alias。close-pair vetoは追加していない。

ROOTにはY閾値・分類定義、`FemtoLambdaMergingDefinition`、EP較正provenance・valid/rejected countersを保存する。EP角×cent9、EPbin×cent9、statusの軽量QAを追加。通常の詳細EP ledgerは無効、candidate storeは明示的に使う場合のみEPD schema2とし、no-EP schema1は維持する。

## 実行入口

|mainconf（config/mainconf/以下）|用途|
|---|---|
|main_auau3p85_anaFemtoLambda_4He_KFParticle_epd8_ylab_merging.yaml|EP8＋YLab分割＋修正merging|
|main_auau3p85_anaFemtoLambda_4He_KFParticle_epd1_ylab_merging.yaml|同じEP valid集合のEP1対照|
|main_auau3p85_anaFemtoLambda_4He_KFParticle_diagnostics_noEp.yaml|EP較正不要のYLab＋修正merging|

再構成は既存`analysis/run_anaFemtoLambda_4He.C`からStChainで行う。mainconf・input/output/jobid・event数を引数で明示する。詳細コマンドは[REFERENCE](../../docs/REFERENCE.md)。既存cent9/Vz/buffer20/両方向mixingは維持した。EP区分はpartner制約に用い、最終raw SE/MEは区分を合算して保存する。SE/MEからのCF・centrality/EP合成法自体を今回変更したものではない。

## 較正runの適用範囲

原run表と較正生成側7表の388runsの全順序が一致する。sourceのcatalog queryは同production/P24iyだがrun≥22158000の部分集合。較正ファイル2本とrun表のSHAを固定した。[生成元・run順・geometry・較正の証跡](../../config/event_plane/xialei_epd_2021P24iy/README.md)。未知runや必要較正cellの欠落は明示errorであり、run番号の置換・近隣runの流用・無較正fallbackは行わない。

全投入snapshot（26994jobs、80982Pico参照、469runs）のうち、73054files＝90.21017%／388runsをカバー。未対応は7928files／81runs。入力listのevent数で重み付けたcoverageは88.27262%で、EP有効event率ではない。24351jobsは全入力が較正対応、2643jobsは未対応runを含む。これは完成26945jobsの部分集合だけの集計ではない。[集計とSHA](../../rootfile/femto_lambda_epd_ylab_merging_validation_20261112_001/epd_run_coverage.json)。全runのEP付き再処理には残81runsの追加較正が必要。YLab/merging修正はno-EP設定で全入力へ使える。

## 検証

SL24y / sl73_x8664_gcc485 / ROOT5.34.38で`singularity_make.sh`を実行し、libStarAnaConfig・libStCommon・libStFemtoMakerをbuild、PASS。新public headerとEpdUtil依存をMakefileへ反映した。全50本の現行Lambda KF mainconfが必須geometry設定を満たす。新hist定義はfemto180、nuclear567。旧key・軸は保持。

- EP回帰20105チェック、YLab/helix回帰462910チェック、cache schema1/2回帰138チェック、全PASS。原Float境界・未知run/SHA/較正key・有効EPD hit・乱数独立性・ME両方向・Y閾値境界・α対He3・磁場/電荷符号・primary fitなし・到達解/長さを確認した。
- 旧fixtureのrun22123040（較正外）はno-EPで4027events、accepted3858、Λ62。既存physics1675hist・39963446cells（flow込み）のcontent/error/Entries/statisticsが完全一致、He4 signal SE9/ME251を維持。Y分割19392cellsの閉包と修正mergingの実ヒスト・alias一致もPASS。
- 較正対応run22158032のactual1992eventsで、EP1/EP8ともaccepted1931・Λ46、EP valid1451/rejected480、8bin全使用、SE853hist完全一致、He4ME63→32。He4SE0だったため次のsampleで追加検証した。catalog記載2229eventsと実Pico1992entriesは別の値として記録する。
- run22159051の先頭10000eventsで、EP1/EP8ともaccepted9761・Λ174。EP valid9374/rejected387、8bin全使用、SE853hist完全一致。実He4 signal SE13で両設定が一致し、ME459→403。He3信号も非空で、MEの282cellsに変化がある。CF形状の物理的改善を判定する統計・検証ではない。
- noEP QA16頁はGhostscript全頁描画と目視確認済み。YLab4質量域の上下SE/MEを追加した。EP QA17頁も出力成功。小統計でnorm域が空の場合はCF未定義と表示する。

[検証study](../../rootfile/femto_lambda_epd_ylab_merging_validation_20261112_001/README.md)、[noEP閉包](../../rootfile/femto_lambda_epd_ylab_merging_validation_20261112_001/noEp_closure.json)、[1万event EP比較](../../rootfile/femto_lambda_epd_ylab_merging_validation_20261112_001/larger_epd_mixing.json)。旧dirty作業ツリーは編集前404filesをbackupし、元からの変更と今回の変更を区別した。過去ROOT・SUMS提出済みsnapshotは変更していない。

## 残Todo

- [x] 今回依頼された解析コード3点の実装・build・小規模実Pico検証。
- [ ] 未較正81runsのrecenter/shift較正を同じ方法で準備・独立検証する。
- [ ] 全統計の再処理と、同じEP有効集合を使ったEP1/EP8の高k*比較。本番投入はまだ行っていない。
- [ ] 元から保留中のmerging用2D ME・符号付き逆rigidity・geometry cache/効率調査。今回のproduction診断修正でこれらを完了扱いしない。
