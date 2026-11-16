# 3.85 GeV Λ–⁴He：新1CEヒストでの高k*系統性

状態：新1CEでの初回系統性走査を完了。高k*低下の原因特定は継続課題。ユーザー指示により過去ログ・Todoを照合し、新しい統合ROOTから比較した。

## 入力とhaddの再確認

2021 P24iyのみ。SUMS `1CE33994716F7896423A80AB994BAB14`、mainconf `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`。
26,945 / 26,994 subjobs、49件除外の検証済み部分統合。2026-11-07 00:40:55 JSTに公開済み。入力2,099,784,951、選択2,066,742,256 events。現在のROOT指紋は公開時と一致した。

11/06 Todoに残った「実行中」は当時の記録。最終判定は `job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14/result.json`、status、archive_resultを優先する。

## 実施する比較

- [x] 過去Todoの³He対照・PID/core-tail・centrality・k*・mT依存と、11/02の未解明baseline問題を確認。
- [x] 新1CEと旧04BFから各1,656個のk* TH1/TH2をROOT5.34/38で抽出。軸・flow・内容・分散・入力不変を検証。
- [x] 信号窓CF、上下半窓、左右sidebandを5統合中央度＋9fine binsで比較。
- [x] nσ窓 ±0.6/1.0/1.4/2.0、正負側を別々に投影。
- [x] rigidity R<2 / R≥2 と1.0/1.5/2.0境界の細区間を別々に投影。
- [x] selectedとprePID、および³He物理質量・³He→⁴He質量仮説を比較。
- [x] mTの0.2 GeV/c²幅区間と規格化域感度、新旧productionの再現性を比較。
- [x] 0.7–0.8、0.8–1.0積分CFと0.5–0.8の傾きを数値化し、図・独立検証を保存。

比較用規格化は0.40–0.50 GeV/c。物理的plateauとして認定しない。raw mass-gate CFを先に比較し、純度補正・LL/FSI fitをこの診断へ混在させない。

## 保存情報による制限

candidate cacheとjoint 8D Sparseは容量対策で停止している。Rとnσの同時条件、低R TOF別CF、任意のmultiplicity/Vz/run/磁場/EP再混合はこのROOTだけでは復元できない。R≥2はTOF必須なのでR区分差をTOF単独の効果と解釈しない。prePIDは既存physical poolを使うQAで、独立した未選別productionではない。

出力先：`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/`。
解析設定・producer・元ROOTを保持し、専用study内に設定・抽出・数値・図を保存する。追加production投入・メール送信・Git commit/pushは行わない。


## 結果

1,680系列（5統合中央度＋9fine bin、3核種表現）を計算。0.7–0.8 GeV/cの積分raw CFは以下。比較用規格化は0.4–0.5。誤差は保存count由来でevent/shared-ME-pair相関は含まない。

|⁴He条件|0–20%|20–80%|
|---|---:|---:|
|全選択|0.90788 ± 0.00328|0.89421 ± 0.00833|
|nσ±1|0.90649 ± 0.00395|0.90111 ± 0.01017|
|R<1|0.85762 ± 0.00515|0.88158 ± 0.01008|
|1≤R<1.5|0.91519 ± 0.00570|0.91281 ± 0.01698|
|1.5≤R<2|0.96986 ± 0.00707|0.98219 ± 0.03329|
|2≤R<5|1.00018 ± 0.01654|1.05448 ± 0.13199|
|2.4≤mT<2.6|0.88991 ± 0.00517|0.89292 ± 0.01067|
|2.8≤mT<3.0|0.96666 ± 0.01057|1.02773 ± 0.05589|

中央試料でnσ±1−inclusiveは−0.00139±0.00222（subset overlapを含む）。own pullの絞りだけで低下は解消しない。rigidity・mT依存が大きいが、R≥2でTOF条件も変わり、Rに相関する混入も除外できない。検出器・TOF起源とは断定しない。

0–20%の³He／³HeAs⁴Heは0.95828／0.95777で、質量変更だけでは⁴Heの全低下を再現しない。Λ上下信号窓は0.90510／0.91145と共通に低下。左右SBは0.96286／0.98188でより平坦。純度補正や新しいBGモデルの採用は行っていない。

旧04BF／新1CEの0–20%は0.90781／0.90788、20–80%は0.89431／0.89421で形の再現性を確認。両productionのイベント重複が大きいため、独立統計の有意差は付けない。0.8–1.0の新inclusive値は0–20%で0.86510±0.00281、20–80%で0.85638±0.00725。規格化域変更だけで傾きは消えない。

## 検証と成果物

新旧ROOTは不変。1,176投影・mass半窓チェックPASS。nσ投影／R投影／既存物理histをflow込みで照合。³He ME overflowに限るTH2F/TH2D丸め差は最大31count、相対2.52e−7。通常binの元の厳しい許容差を維持し、flowだけFloat32 ULPで別検査。0≤k*<1のCFへ影響なし。

別実装で直接native積分した318数値・subset共分散と24有限差分導関数検査がPASS。描画84,000binのうち欠測7,804はNAのまま表示し、ゼロに置換しない。低統計sliceの傾きには実際に使用したbin範囲を保存。

- [図と結果の報告](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/lambda_he4_highk_systematics_20261111.pdf)
- [18頁比較atlas](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/cf_highk_systematics_20261111.pdf)
- [条件定義・解釈・再現手順](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/README.md)
- [数値表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/metrics.tsv)
- [主検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/validation.json)／[独立検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_highk_systematics_20261111_001/independent_validation.json)

## 残る解析

次に切り分けるのはR/mT compositionを合わせた比較と、候補/PID固定のcontrolled remix。³He core→⁴He窓の混入率定量、最終純度補正・LL/FSI fit、mixing条件差の原因検証は未完了。candidate記録なしのhist-only試料だけで解決したとは扱わない。今回、追加production・設定変更・メール・Git commit/pushは行っていない。

最終成果物は説明2頁＋atlas18頁の計20頁。数値・ソース指紋、PDF頁数・EOF・組版、代表頁の目視、新旧入力ROOTの最終指紋不変を確認し、studyの`completed.json`へ記録した。


## CF前のSE/ME k*重ね描き：作図進行 — 09:37 JST / 2026-11-12

ユーザー依頼により、CFを計算する前のselected inclusive Λ–⁴He信号SEとMEを同じpadへ重ね描きする。主図は0–20%、対照として20–80%と0–80%を作成中。元のSE countは変更せず、MEだけ一定係数SEnorm/MEnormでscaleし、0.40≤k*<0.50 GeV/cの積分を一致させる。比較用の表示規格化であり、各binの一致や物理的plateauを仮定しない。

入力は上記と同じ2021 P24iy・SUMS `1CE33994716F7896423A80AB994BAB14`、26945/26994成功の検証済み統合ROOT。mainconfは`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`、configlogは`job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14.txt`。物理TH2 `hKstarSEVsCent_lambda_he4_signal` / `hKstarMEVsCent_lambda_he4_signal`を中央度投影する。Λ信号窓1.10708–1.12484 GeV/c²、追加⁴He PID区分なし。背景差引き・純度補正なし。

native k*幅0.005 GeV/cを4binまとめ、表示幅0.020 GeV/cとする。規格化のnative積分と予定ME係数は以下。

|中央度|SEnorm [0.4,0.5)|MEnorm [0.4,0.5)|MEscale=SEnorm/MEnorm|
|---|---:|---:|---:|
|0–20%|122358|4461431|0.027425729547313406|
|20–80%|18213|666345|0.0273326880219706|
|0–80%|140571|5127776|0.027413638973309288|

保存先は`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/`。この段落は09:37時点の進行記録。作図・主検証・保存ROOT再読の完了は後の作図完了節を参照。

- [x] 主図`se_me_cent0_20.png` / `.pdf`、対照`se_me_cent20_80.png` / `.pdf`と`se_me_cent0_80.png` / `.pdf`。
- [x] 3頁の`se_me_overlay.pdf`、`overlay_histograms.root`。
- [x] `settings.json`、`normalization.tsv`、`validation.json`、`README.md`。
- [x] `independent_validation.json`、`completed.json`に最終検証・作図完了の証跡を記録。
- [x] [study記録](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/README.md)へ入力不変・投影・rebin・規格化・PDF頁数と目視確認を保存し、後の節へ作図完了を追記した。

追加production/QA job・本番source/config・mailは変更しない。track-mergingの追加Pico再読はユーザー指示で保留中。完了済みの高k*系統性比較やtrack-merging最終結果は保持する。


## CF前のSE/ME k*重ね描き：作図完了 — 09:50 JST / 2026-11-12

3中央度のPNG/PDFと保存ROOTを作成し、主検証と保存ROOT再読まで完了。主図は0–20%で、各図左を0–1、右を0.3–1 GeV/cとした。SEは元countの黒点、MEは前節の一定係数を掛けた赤線、灰色帯は規格化域0.40≤k*<0.50 GeV/c。native .005をrebin4して表示幅.020 GeV/c、smoothingなし。背景差引き・純度補正なし。

- [主図0–20% PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_cent0_20.png)／[PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_cent0_20.pdf)。
- [20–80% PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_cent20_80.png)／[PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_cent20_80.pdf)。
- [0–80% PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_cent0_80.png)／[PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_cent0_80.pdf)。
- [3頁PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_overlay.pdf)／[保存ROOT](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/overlay_histograms.root)／[条件と再現手順](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/README.md)。
- [settings](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/settings.json)／[規格化表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/normalization.tsv)／[主検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/validation.json)／[独立検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/independent_validation.json)／[完了記録](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/completed.json)。

ROOT5.34/38の`SEME_SUCCESS pages=3`を確認。physical TH2と前回CF用QA投影の通常200binはSE/MEとも全てexact一致。規格化・一定scale・rebin count/errorのclosureは全中央度PASS。元1CEのsize/mtime/inode/deviceとUUIDは不変。保存ROOTはTH1D18・TParameter3・TCanvas3の24objectを再読し、`OUTPUT_ROOT_REOPEN_PASS keys=24`を確認した。

PNG3枚は目視良好。結合PDF3頁と個別PDF各1頁をGhostscript nullpageとPage object数で検証した。初回runnerはROOT作図成功後にpdfinfo不在で後処理だけ停止したため、インストール済みGhostscriptへPDF確認を切り替えて復旧。ROOT再描画・図や数値の変更は不要だった。独立TSV全bin検証もPASS。主検証・独立検証・保存ROOT再読により作図の完了を確認した。

規格化域では積分が一致し、約0.6 GeV/cから高k*側ではscaled MEがSEより上になる。過去の高k*CF低下と整合する分布形状差であり、この図だけでSE/MEどちらが歪んだか、mixing・PID・mergingのどれが原因かは確定しない。scale由来の不確かさやevent/共有ME pair相関を追加評価したものではない。追加production/QA job・本番source/config・mail変更なし。track-merging追加2D ME/逆rigidity/geometry再読の保留を維持する。


## ME作成方法の本番照合 — 10:10 JST / 2026-11-12

保存runtimeと現行実装の一致を確認。MEはcent9×Vz別の直近最大20候補付きeventを両方向総当たり。0–20%を一括poolにはせず、固定event間隔なし、pool1件から開始、1/pool補正なし。d/tだけのeventもpool枠に入り、job内run/file境界reset・run/B/EP別poolなし。Vz198–202cm/10区分はeventcutでなく範囲外clamp。実際のevent差・端pool占有率は未評価。今回source/config/ROOT変更・job投入・mailなし、track-merging追加確認は保留。

- [x] 1CE投入設定・frozen sourceに基づく[ME手順の監査](../../../../analysisnote/20261112/lambda_he4_me_construction.md)。
- [ ] 端Vz pool占有、実event gap、有効Λ/⁴He partner数、pool立上がり・event重みを診断し、event gapとbuffer深さを分けたcontrolled remixを比較する。まだ再投入していない。


## φp/φdのMEとの比較 — 17:49 JST / 2026-11-12

- [x] [StPhiMakerとφtree下流のMEを比較](../../../../analysisnote/20261112/phi_p_d_me_comparison.md)。φp/d CFは同cent9×Vz×EP、全accepted eventがpoolへ入り、Vzを実cutする。固定gapなし・満杯待ちなし・1pair1countはΛと共通。StPhiMakerのKK質量背景MEとは別処理。
- [ ] φ最終runの実効buffer/mode/Vz値は保存runtime/argv/ROOT metadataで確定する。現在は標準参照YAML・φ提出記録が欠落。残るrndsamp設定とwrapper既定値を実際の採用値とは区別する。


## oura／xjiang2共有sourceとの比較 — 22:13 JST / 2026-11-12

- [x] ouraの提出済みφ設定・9/19/9/28 job argv・保存MEsourceを確認。標準Vz実cut/10×cent9/EP1/buffer20/bufferAll。全accepted pool。
- [x] Xialei Jiangのaccount xjiang2を確認し、dΛ系とφp保存codeを探索・比較。dΛ系はEP8/Vz8/depth50/random replacement/単方向。φp保存codeの同event除外等の注意点は最終採用版か未確認と区別。
- [ ] Vz採用・pool入場・EP・更新法/深さ・方向を独立に変えるcontrolled remix。今回job投入なし、merging追加QAは保留維持。

[比較報告と制限](../../../../analysisnote/20261112/external_me_comparison_oura_xjiang2.md)／[source/config/argv指紋](../../../../analysisnote/20261112/external_me_comparison_evidence.json)。


## EP導入と分割CF合成の検討 — 22:33 JST / 2026-11-12

- [x] 現macroのcount合算→group規格化→比を確認。Xialei EPは混合制約に使い、EP合算rawhistへ保存することも確認。
- [x] 区分別ME規格化後の合成との厳密な重み対応、固定CF平均との違い、誤差/covarianceを整理。既存1CEの合成方式差はhigh-k約0.0005–0.0015で低下を解消しない。
- [ ] EP角のharmonic/較正/自己相関/invalid集合を定義し、EPなしと粗い/細かい区分を同一event集合で比較。既存補助TPC Q2は未較正、1CEの後付け分割不可。
- [ ] cent9×EPのSE/MEと有効mix量を保持した再処理・job/block covariance評価。本検討では新job投入なし。

[検討報告](../../../../analysisnote/20261112/ep_and_cf_combination_review.md)。LambdaLegacyはSE/ME保存のみでproducer保存CFなし、今回CF規格化は後処理0.4–0.5という点も同報告に明記した。

合成方式比較の[独立source・設定・全数値・検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_centrality_combination_20261112_001/README.md)を保存済み、全PASS。規格化域SEペア数の固定重みCF平均も別の観測量として比較し、追加global再規格化を加えていない。


## Xialei方式EPの採用意向・定義確定 — 23:13 JST / 2026-11-12

- [x] 上流生成元・2021run表・較正file所在を確認。EPD East AB(row1–8)の較正済みΨ1、8分割を採用対象とする。TPCΨ2とは別。
- [ ] 較正ROOTのkey/axes/entriesと対象run coverage、較正用cent・元valid条件を照合。
- [ ] EpdHit読込、EPD provider、store/source/schema/guardとbin定義を移植し、元codeとのpilot照合を実施。現時点では未実装・新jobなし。

[採用仕様と残検証](../../../../analysisnote/20261112/xialei_epd_method_adoption.md)。


## EP・YLab・merging解析コード反映 — 01:49 JST / 2026-11-13

- [x] Xialei EPD一次EP8と同有効event集合のEP1対照を実装。実Pico10000eventsでSE853hist完全一致、実He4SE13・ME459→403、8bin全使用を確認。
- [x] 物理He4のYLab strict −0.5分割でSE/ME32histを追加し、両ME方向・mass域・inclusive閉包を検証。
- [x] production track-merging診断をSTAR helix/実B/global rigidityへ修正し、既存物理SE/MEの不変を確認。
- [ ] 残81runsのEP較正を準備・検証し、全統計再処理へ進む。現在の較正は全投入入力filesの90.21017%／388runsだけを覆い、未対応runへ流用しない。

[実装入口・検証・残Todo](../../../../analysisnote/20261113/lambda_he4_epd_ylab_merging_implementation.md)。本番SUMS再投入・mailなし。過去のCF合成方式やclose-pair vetoを変更していない。元のmerging追加調査の保留は維持。


## EPの全run対応方針の補足 — 10:49 JST / 2026-11-13

[既存係数の共有可否・Pico内情報・EP確認図](../../../../analysisnote/20261113/ep_calibration_transfer_and_all_runs.md)。未対応81runは前期day121–125。今のrun別EP表に係数が無いという状態であり、81runすべて個別に新規較正する以外の方法が無いと断定しない。

- [ ] 対象81runで明示した共有係数候補を独立標本で検証する。
- [ ] 共有が適合しなければ前期の期間共通／run別recenterと、それに対応するshiftを準備・独立検証する。
- [ ] 全runのEP1/EP8再処理。今回source/configの有効化・SUMS投入は未実施。
- [x] 補正後EP角・実8区分占有率を確認する境界線入りPNG/PDFを作成（1万event pilot、有効9374）。


## ユーザー方針更新 — 18:56 JST / 2026-11-13

未較正81runの解析・共有補正試験・新較正作成は保留とする。まず較正済388runでEP分割の有用性を確認し、有用と確認できたら81runの較正を行う方針。前節の共有補正検証／追加較正案は再開待ちであり、現在の実行予定ではない。
