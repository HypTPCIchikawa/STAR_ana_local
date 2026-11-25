# Λ–He Δφ*–Δη作成・分布の監査

状態：監査完了。作成側の不具合を確認。本番ソース・設定・ROOTは未変更。修正実装・再生成・cut採用は未実施。

ユーザー依頼：他解析と異なり原点ピークが見えないため、実ヒストと作成ソースを確認。

## 入力・再現性

2021 P24iy、SUMS `1CE33994716F7896423A80AB994BAB14`、26,945/26,994の検証済み部分統合（49件除外）。mainconf `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`。runmetaと同IDのconfiglogを照合。

提出時runtime archiveと現行の`FemtoLambdaLegacy.cxx/.h`、StFemtoMaker、Maker YAML、核hist YAMLの5ファイルがbyte一致。解析側だけの後変更ではなく、投入時から存在する問題と確認した。

## 確定した問題

- **曲率補正の符号**：`FemtoLambdaLegacy.cxx:997–998`のEventField版も、φp−φHeへ+asin(Ap)−asin(AHe)を加えている。STAR helixでの位置角はφ(r)=φ0−asin(c qBr/(2pT))。±B・±電荷・複数pT/φ/pzの192ケースで確認。正確な係数の−asinは最大誤差1.33e−15 rad、現行+asinは0.86657 rad。真に同じ位置の理想ペアが現行式で+0.17935radへ動くため、横軸全体の正負規約差ではない。
- **Λ娘へのprimary fit要求**：同:987–988がpPt()>0とpMomを使う。primary情報がないsecondary娘を落とし、CFとは違うsubsetになる。SL24y StPicoTrack.h:72/178で定義を確認。同解析の娘acceptanceは:848–850でgMomを使っている。
- **旧名ヒストのHe曲率**：同:992–993は固定B=+0.5T、raw rigidity pTにZ=2を掛けるためHeの曲率項が過大。物理charge2/pTphys=2/(2pTraw)でZは相殺する。EventField版はこの因子とイベント磁場を改善するが、上2点が残る。
- **半径到達性未検査**：r=1.4m固定、asin引数をclampするだけ。FXT Vz≈200cmと娘originを考慮して、両trackが共通半径へTPC内で到達するかを調べる必要がある。

trackIDを配列indexに使う旧問題は修正済み。現存の6histは正しい配列indexを参照する。

## 実ヒストで見えたこと

6個＝⁴He/³He/³HeAs⁴He ×旧名/EventField。Λ信号SE、proton娘のみ、全k*・全centrality積分。ME、π娘、k*/R/cent角度軸はない。両軸200bins[-0.25,0.25]、幅0.0025。

⁴Heの全信号SEは2,641,987（k* flow含む）。QA Entriesは両版2,253,813＝85.3075%。表示領域内は旧58,100／EventField57,813で約2.6%、平均約1.45count/bin。約97.4%が角度flowに入るのは狭い範囲の性質で、データ破損・削除ではない。また、fill率85.3%はQA admissionで、検出器が14.7%のpairを失ったという意味ではない。

|原点box|旧名|EventField|
|---|---:|---:|
|±0.005|19|18|
|±0.01|89|80|
|±0.02|368|320|

拡大・投影でも鋭い原点過剰はない。EventField原点±0.01／隣接square-ring密度比1.00±0.13（count統計のみ）。誤った角度定義のraw SEなので、真のmerging効率を判定したものではない。

全6個でEntries=全bin積分。旧名/EventFieldのEntries一致、³Heと³HeAs⁴Heの軸・内容・誤差・Entries完全一致を確認した。

## CF・物理解釈

mergingは近接ペアの再構成効率低下であり、raw SE原点ピークを必ず生じるものではない。splittingと区別し、SE/MEやpair効率で見る。[ALICE一次資料 §2/§3.3](https://arxiv.org/pdf/1507.06842)。他実験のcut値は採用しない。

物理CFは同ソース901–921行で先にfillされ、mergingは926行のQA。現設定はclosePairEnabled:false。このQA誤式が既存CFに直接cutを掛けているわけではない。実際のpair効率が高k*低下に寄与するかは未評価。

## 成果物・残作業

- [11頁監査レポート（説明2頁＋図9頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_audit_20261111_001/track_merging_audit_20261111.pdf)
- [日本語の詳細・再現手順](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_audit_20261111_001/README.md)
- [STAR helix検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_audit_20261111_001/helix_sign_validation.json)／[実ヒスト数値](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_audit_20261111_001/hist_audit.json)

修正方針はglobal track幾何、正しい曲率符号・charge/rigidity・イベント磁場、共通TPC内到達半径。proton/π娘両方、SE/ME両方、k*/cent依存も保存する。STAR helixと小統計PicoDstで選択を固定した検証後に、本番再生成・close-pair cutを検討する。

今の2Dヒストを定数shiftや全体軸反転して修復できず、QAから除外した娘も復元できない。candidate cacheにも必要なglobal daughtertrack幾何がなく、元PicoDstのtrackへ戻る必要がある。今回、修正実装・追加SUMS投入・メール送信・Git commit/pushは行っていない。


## 続報：全再投入を避けるQA専用再計算（2026-11-11 JST）

ユーザーの「このヒストだけを短時間で見たい」に対し、保存候補＋対応PicoDstのEvent/Trackだけを使う独立StudyMakerを作成。KF・PID・イベント選別の再実行なし、本番source/config/libraryは変更しない。

ローカル4027-event Picoと3858accepted候補cacheで成功。候補が必要な52イベントだけPicoをrandom readし、StChain処理＋入力読取1.83457秒（コンパイル・ROOT起動時間除く）。全旧6histの200×200binとflow、選択済みSEペア数の244,827項目が完全一致。新しいglobal負asin、actualhelix/TPC内1.4m、共通到達半径の走査版をproton/π娘両方に作成した。

[QA専用ROOT・96hist](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_quickcheck_20261111_001/smoke/qa.root)／[実行・入力・ソース指紋検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_quickcheck_20261111_001/finish_validation.json)／[実カウント](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_quickcheck_20261111_001/smoke/inspection.json)。

本選択He4は9pair、He3は29pairでピーク有無を判定できる統計ではない。He4 GlobalMinusAsinは9pair/表示範囲1/原点±.0125は0。固定1.4m実helixはproton7/π9pair。共通半径で両娘9pair。TPC PID前のgeometry controlは3044pair（負電荷325を含む別population）、proton固定helix2299/表示範囲70/原点±.0125は2。controlの形を最終He4 sampleと混同しない。

既存2AF17000にはHe4SE115/He3SE474、26993には101/493があり、正常候補ROOTと元3Picoの順序付きlistが保存されている。全jobを回す必要はなく、この方式を少数jobへ拡張可能。ただしROOT5直接remote入力はNetx不足で使えず、同版pluginの隔離loadもXrdUtils.so.2不足。17000の元3Pico計2.62GBの追加staging・再計算は今回行っていない。今後は元Picoをローカルに用意できる標本を使う。smokeの1.8秒をremote入力・大標本へそのまま外挿しない。

angular SEのみなので効率損失や物理CFへの影響は未判断。近接pair人工試験の原点再現と、実データにピークがあるかは別の検証。新規描画sourceはauto-reviewがAGENTSの明示承認規則を理由に拒否したため未作成。既存ROOT描画によるderived画像出力に切り替え、本番処理は変更していない。


短時間テストの確認図は、拒否された新規描画sourceを作らずROOT標準のTFile/TH2::Draw/TCanvas::Printで生成できた。
[修正版の確認図PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_quickcheck_20261111_001/smoke/quick_merging_overview_v2.pdf)／[PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_quickcheck_20261111_001/smoke/quick_merging_overview_v2.png)。selected-proton近似式、selected-proton実helix、selected-π実helix、全電荷prePID control実helixを固定r=1.4mで表示。新しい物理結論は付けず、計算・読出しの動作確認として提示する。

最終確認図v2は全panelのcount色尺度0–2を共通にし、空histの負の自動色尺度を解消。図を目視し、完了指紋をquickcheck studyのcompleted.jsonへ保存した。


## 優先Todo：本番解析コードのデバッグとQA統計の拡張

ユーザーの明示指示（2026-11-11 JST）：可能な限り元PicoDstを読み、QA専用の再計算で統計を増やす。また、**解析コードそのもののデバッグを行うことをTodoへ明記する**。独立StudyMakerの動作確認と、本番解析コードの修正完了を区別する。

- [x] **本番解析コードのデバッグ**：`StMaker/common/FemtoLambdaLegacy.cxx::FillMerging`の曲率符号、Λ娘のpMom/pPt要求、He電荷/rigidity、実磁場、secondary origin、TPC内の共通到達半径を修正・検証した（2026-11-13 JST、[実装・検証記録](../../../../analysisnote/20261113/lambda_he4_epd_ylab_merging_implementation.md)）。全統計再投入・merging効率確定は未実施。
- [ ] 修正前の旧6hist・物理SE/ME・選択候補が再現されることを小統計で確認し、QA修正によって物理ペア選別やCFが変わらないことを確認する。
- [ ] proton/π娘の両方、SE/ME双方、k*/centrality依存を調べられるQAの設計と検証を行う。close-pair cutの採用は、その後に別途判断する。
- [x] 保存候補＋元PicoDst Event/Trackだけを読む独立QA処理を小統計で検証（52イベント読取、244827一致検査PASS）。
- [x] 2021 P24iy候補ROOTと元PicoDstの順序付きlistを照合し、利用可能標本の一覧・重複・対象coverageを確定（primary 10,216＋ROOT完走検証済みsecondary 9＝10,225 jobs、30,675 unique PicoDst）。最終採用には各workerの完全検証を要求。
- [x] QA専用runtime/configを固定し、worker-local scratchへ入力を一時stageするfarm pilotを実行。正しいイベント対応・完走・旧hist一致・読取時間・容量を検証する。
- [x] pilot結果と利用可能資源を見て、制限した同時実行数で利用可能標本へ拡張する。原Pico・候補ROOTは保持し、削除はworkerが自ら作成した一時stageコピーに限定する。
- [x] 成功分だけを重複なく集計し、実測イベント数・本選択ペア数・原点count・図・除外/失敗一覧を更新する。

統計拡張の出力先：`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/`。新しいfull KF/PID解析や別datasetへの変更ではなく、元の選択を固定したmerging QAの読取りを進める。pilot完了や全対象完了は実測前に主張しない。


## 統計拡張の投入・継続実行 — 14:56 JST / 2026-11-11

利用可能な保存候補全体から、**10,225 jobs / 30,675 PicoDst**を対象に確定した。全て2021 P24iy、URL・LFN重複0。対象の元producer実測inputは823,337,640 events、accepted candidate framesは812,995,885。これは対象coverageであり、今回QAで読み終えたevent数ではない。primary 10,216件に、最終stdoutは欠けるがROOTのcompleted・非Recovered・UUID・provenance・tree/counter一致を検証した9件を加えた。各workerでpayload全読取、元Pico全Entries、対応event/run/PV/B/track、旧QA全bin/flow・SEペア数の再検証を通過した出力だけ採用する。

- [対象inventory](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/dispatch_inventory.json)／[独立照合](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/dispatch_inventory_validation.json)
- 固定runtime v2はSL24y/ROOT5.34/38でbuild・local smoke PASS。local fixtureのHe4 9 pairはfarm増加分に数えない。
- 初回pilot `708262` はXRootD serverがchecksum queryを提供せず転送開始時点で失敗し、新規QA増加0。失敗出力・初回runtimeは保持。通常xrdcp＋前後のremote size/MTime、local SHA256/adler32、ROOT内容の完全検証へ変更したv2を別保存した。server checksum確認済みとは主張しない。
- v2の実データpilotは **Condor `708273`**、元subjob 17000/26993。現時点は両方Idle、開始0・完了0。要求条件の不一致は見つからず、CAS group割当/優先度による待ちと見られる。
- 両pilotの厳密検証に合格した後に、準備済み10,223 jobsを一度だけ自動投入する。`max_materialize=32`、`max_idle=16`、worker-local scratchへstageする。現在は全量factory未投入。
- tmux `he4_merging_dispatch_20261111` がpilot合格待ち、`he4_merging_campaign_20261111` が成功receipt検証・重複なしincremental hadd・原点count集計を継続する。lockとsubmit intentで二重投入を防止。原入力は保持し、一時stageコピーだけworkerが削除する。

[最新の実処理状況](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/status.json)／[投入gate状態](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/advance_status.json)。成功分のhaddと数値は同studyの`aggregate/latest.json`に保存される。source input coverage、cache frames、実Pico Event/Track読取数、selected He4/He3 pair数を区別して更新する。

**本番解析コードのデバッグは未完了**。独立QA処理の検証成功と混同しない。統計拡張後に修正QAの形を比較し、本番FillMergingの修正・小統計検証・SE/ME/両娘QA整備へ進める。


## ローカル先行標本の追加 — 15:17 JST / 2026-11-11

farmのfairshare待ちと並行して、保存済みPicoDstを探索し、既知4027-event sample以外に正式名の3本を確認した。0BEF/04BFの保存SUMS listと元URL、remote size、ROOT nonZombie/nonRecovered・fEND・PicoDst/Event/Track/TrackCovMatrixを照合済み。2021 P24iy、実測Entries34,603＋34,841＋13,140＝**82,584**（catalog98,606とは別）。既知sampleを含む4LFNはfarm対象30,675LFNと重複0。

対応する2AF候補は保持されていないため、この3本だけを既存candidateStoreWriteNoJointSparse設定で個別・1CPU逐次処理する。過去の元jobは別の不完全companionも含み未完走だったが、正式3本のpayload健全性は今回の完走で別途確認する。producerライブラリ14本は、11/02のbatch-matched成功smokeのSHAと全一致。nice10・1本20分上限で元入力や本番source/config/libraryを変更せずに実行中。

[入力来歴・実測Entries](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/input_provenance.json)／[先行読取りの進捗](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/producer_status.json)。各producer完走後、別に固定したlocal用QA runtimeで全bin一致・イベント対応・入力不変を検証する。成功分を既知標本と合わせて図を作り、farm実績とは別に記録する。現時点では先行読取り中で、QA完了や増加ペア数は未確定。


## ローカル先行標本の結果（完了）— 15:33 JST / 2026-11-11

追加3Picoの82,584eventsは全てproducer完走・元入力不変・既存lib不変で成功。既知4027-event標本も同じlocal QA runtimeで再検証し、全4 sourceで244,827項目ずつの旧hist全bin/flow・SEペア数一致、UUID/completed/nonRecovered、元Pico Entries=eventLedger=inputEvents、accepted cache数・選択Λ数を確認した。config13ファイルの本番設定と固定snapshotのSHA一致も確認。

| 区分 | source input coverage | accepted cache frames | QAで実読取したPico events | He4信号SE | He3信号SE |
|---|---:|---:|---:|---:|---:|
| 既知標本 | 4,027 | 3,858 | 21 | 9 | 29 |
| 今回の追加分 | 82,584 | 79,086 | 447 | 115 | 502 |
| 合計 | **86,611** | **82,944** | **468** | **124** | **531** |

今回、producerは追加82,584eventsを全読取した。QAは候補が必要な468eventだけEvent/Trackを読取り、82,944cache frameは全decodeした。coverageとQA実読取数を混同しない。He4は9→124で約13.78倍。farm実績への加算は0、main対象と4LFNの重複0、除外0。

- [x] local追加3本の全イベント読取り・候補生成。
- [x] 全4 sourceの厳密QA・入力対応・旧hist closure。
- [x] 成功4本のhadd、全96hist/flowを含むcontent・Entriesの正確合計一致、error²を含む1,699,968値の検証。Sumw2未保存のerror²はsqrtの丸めのみ最大8ULP許容（758値、最大差1.42e-14）。整数countの検査は緩めていない。
- [x] 修正版proton/π娘と旧EventField比較の4頁PDF・4PNGを生成、入力ROOT不変、PDF/図を目視確認。

[4頁確認図PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/plots/jobs_00004_e46f8399caba_ba5ba010af/merging_maps.pdf)／[最新結果と実統計](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/status.json)／[厳密QA検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/runtime_validation.json)／[haddと数値](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/aggregate/latest.json)／[図の検証・画像一覧](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/plots/latest.json)。

He4 GlobalMinusAsinの表示範囲±0.25内はproton5／π7pairで、全methodの原点±0.0125および±0.025は0。**原点ピークの有無を判断できる統計ではない。** raw SEでありMEとの比較・効率評価は未実施、ClosestSpatialEtaはpair毎の最接近半径という条件付けも異なる。今回の成功は計算手法・入力対応の検証と統計増加であり、本番FillMergingデバッグやclose-pair cut採用の完了ではない。

大規模対象10,225jobs/30,675Picoのpilot708273は依然Idle。両pilot合格後の10,223-job自動投入、成功QAの自動hadd、図の更新はtmuxで継続中。本番解析コードのデバッグは上記優先Todoの未完了項目として維持する。

独立の最終検査もPASS：[元入力・候補/QA ROOT・lib/config・全数値の再照合](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/local_extension/independent_results_validation.json)。最終図はsource Pico entriesとQA実読取数を明記した版を採用した。


## 統計倍率の見積りとfarm転送条件の更新 — 15:58 JST / 2026-11-11

ユーザーの30,675ファイルで何倍になるかという質問に対し、元producer実測823,337,640events / local標本86,611events＝**9,506倍**、accepted framesでは9,802倍と確認。124 He4 SE pairsからの単純外挿は約118万。実際に完了したfarm140jobsは11,512,341events / He4 14,627pairsで、これを用いると全対象で約105万pairs（124の約8,436倍）。**概ね100〜120万pairs、約1万倍**という見込みであり、検査通過分と入力分布に依存する。原点付近のcountは現local0countから直接外挿せず、実ヒストで確認する。[算出根拠](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/statistics_projection.json)。

v2 pilot2本が合格してmain factory708357が開始した。ところがredirectorが同じURL/sizeでもIdと過去MTimeの異なる複製を返し、更新時刻の厳密一致で転送を拒否した。停止操作でJobMaterializeLimit=0を停止と取り違え一時展開が増えたため、JobMaterializePaused=1とLimit32へ訂正。factory Held/Next568、実queue0を確認し、旧出力とitemdataを保管した。

568 main＋2pilotの最終内訳：成功140、MTime判定失敗400、厳密な入力coverage/identity失敗23、source/local size不一致7。後者30件は検査を緩めず除外。MTime400＋未展開9655＝10055件をimmutable runtime_v3へ割当てる。v3はMTime/Idを監査情報とし、before/after/local size、digest、ROOT fullcoverage/event identity/旧QA closureを維持。physics source/include/config/libはv2と全SHA一致。実STAR QA付きの回帰検査・不正size拒否とcleanupもPASS。[v3検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/runtime_validation_v3.json)。

v3の実pilot708368（15642/15660）も両方completed/244827 closure PASS、He4合計139pair。残り10053件の再投入を準備する。v2成功140件は保持し、同じsourceを重複集計せず、source別保存先/runtimeSHA/policySHAを照合する監視へ切替済み。[source別ルート検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/result_routes_v3_validation.json)。local標本とfarmは別集計、4LFN非重複、本番コードdebugは未完了のまま。


## v3 main再開・実factory確認 — 16:16 JST / 2026-11-11

独立47項目のcohort検証と両v3 pilotの厳密QAを通過し、rootのready-file gateが07:15:01 UTCに**main cluster 708399 / 10,053 jobs**を一度だけ投入した。400件のMTime再試行＋未展開9,655件から成功pilot2件を除いた集合で、元入力行・URL・全hashを照合済み。旧factory708357の全itemdata、Proc0..567のhistoryと未展開全IDは削除前snapshotに保存した。

07:16:06 UTCの実factoryで**JobMaterializeLimit=32 / MaxIdle=16、PRESENT32・RUN32・IDLE0**を確認。各実jobの1CPU/2GB/8GB、runtime_v3/policy、予定subjob所属もPASS。v2成功140＋v3 pilot成功2の計142sourceは保持し、厳密なcoverage/identity・size不一致30sourceは除外を維持する。source別runtime/policy/store検証により重複採用せず、成功QAの自動haddと確認図更新を継続する。local標本とfarmは別集計、本番コードのdebug完了とはしない。

[実投入receipt](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/main_v3_submission.json)／[実factory・resource検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/main_v3_launch_validation.json)／[再投入cohort独立検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/recovery_v3_inventory_validation.json)。


## 完了見込みとメール通知 — 16:32 JST / 2026-11-11

ユーザーが完了見込みと終了後メール通知を依頼。v3 main708399の実測のみでETAを算定：07:23:54 UTC時点の直近2/3/5分は23.5/24.0/25.4 jobs/min、成功worker wall平均64.9s、残り9849。32枠が継続し転送速度が同程度なら、job終了は日本時間11/11 23:19中心、最終hadd・検証・図に仮15分を加えて**23時半前後、幅を見て23時〜11/12 1時ごろ**。資源割当や転送時間の変動で前後する。[ETA算定根拠](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/eta_estimate_20261111.json)。

- [x] 完了メールを`yudai.ichikawa.d3@tohoku.ac.jp`宛てに設定。userの今回の明示依頼に基づく。現在は未送信。
- [x] tmux `he4_merging_mail_20261111`で常駐、120秒間隔。全sourceの現在marker/receipt/summaryを再読、成功/失敗集合の全件確定・queue空・最終aggregateの集合/counter/ROOT SHA・それに一致するPDF/PNGのbytes/SHAを確認してから一度だけ送信。
- [x] 既存`script/watch_job_notify.py::send_once`を再利用。MTA呼出し前にattemptを永続化し、送信結果不明時は自動再送せず二重送信を防ぐ。失敗・除外があれば件名と本文に明示。成功0件でも結果なしとして通知可能。
- [x] 通知fixture/mock MTAの14検査PASS（実メールなし）、現在の未完了campaignでは送らないことを確認。plot monitorの最終snapshot競合も修正・15検査PASS、常駐を再起動済み。

[通知設定](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/notification/settings.json)／[常駐と検証記録](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/notification/launch.json)／[通知の現在状態](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/notification/watch_status.json)。メールには成功/失敗件数・He4/He3ペア数・最終PDFと結果保存先を含める。メール配送はlocal sendmailの受付記録を残し、受信箱への到達と同一視しない。


## 11/12未明の進捗とETA更新 — 00:30 JST / 2026-11-12

userの完了確認に対し、00:26 JSTでは未完了。検証済みstatusは成功9177/除外127/残921、He4信号SE922798、成功9177件のhadd/図が同じSHAで更新済み。statusは6〜7分のpoll時間差があるためfresh Condorで再計算し、00:27は未展開802＋queue32＝残834、31Running/1Idle/Held0を確認した。

直近10/30分は15.7/16.27jobs/minで、初期24〜25/minより低下。32枠が続けばjob終了01:15〜01:40、最終hadd/図/メールは11/12 **01:30〜02:00 JST**が更新目安。前の23時〜翌1時見積りを更新する。campaign/plot/mailの3daemonは稼働、integrity/merge error無し、mailはcampaign_not_terminalで待機し送信attempt無し。全job終了と最終ROOT/図検証後に1通送る設定を保持。今回はsource/config/runtime/queue変更・取消・再投入・メール送信なし。

[ETA根拠](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/eta_update_20261112.json)／[進捗snapshot](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/progress_20261112_002.json)。


## 現時点統計のΔφ*–Δηテスト図 — 00:37 JST / 2026-11-12

userの指示で全job終了を待たず、検証済み9272sourceの固定snapshotを使い確認図を作成。元source入力741,914,173events、He4–Λ信号SE **932,826pairs**。ROOT SHA `0712dd157070f18ab5dc1dfe69c9bc9a05ab16f021dd3d0523562af4367bd12d`、途中統計で最終完成と混同しない。

- [固定r=1.4mの陽子/π、全域＋原点zoom、1頁PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/quicklook_20261112_003/focused_fixed_radius/merging_fixed_radius.pdf)／[PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/quicklook_20261112_003/focused_fixed_radius/merging_fixed_radius.png)
- [旧EventFieldと複数修正定義を比較する4頁PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/quicklook_20261112_003/jobs_09272_0712dd157070_ba5ba010af/merging_maps.pdf)
- [入力・作成検証と原点数値](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/quicklook_20261112_003/completed.json)

固定1.4mのSTAR helix＋global momentum ηでは、protonのQA採用858058/表示±.25内23863/原点両軸±.0125は26、πは903265/26000/71。全93万pairと狭い表示域/原点countを区別する。原点squareと隣接2倍square ringの面積補正密度比は、halfwidth.0125でproton0.79/π1.09、.025でproton0.56/π1.02。固定半径とGlobalMinusAsinに鋭い原点過剰は見えず、proton中心密度は低め、πは概ね周囲と同程度。rawSEの形の記述で、MEなしの効率評価やcut採用の結論にはしない。ClosestSpatialEtaはpair毎の最接近半径という条件付けがあり、ピーク根拠に使わない。

source ROOT/hash・全bin/count一致など既存検証と図出力後の入力不変、PDF4頁＋1頁・PNG・目視を確認。全域は表示用40bins、zoomはnative80bins、smoothingなし。1頁図は各panel独立color scale。raw信号SE、全保存k*/centrality、background subtractionなし、新close-pair cutなし。farm/自動hadd/最終図更新/完了mail通知、本番解析コードは変更していない。


## 全job終了・全成功統計の最終Δφ*–Δη — 07:15 JST / 2026-11-12

userの全統計図作成依頼に対応。全10225sourceは終了、**10087成功＋138除外、pending0**、Condor708399のqueue/factory空を確認。成功分の最終haddは11/12 01:29 JST、最終比較図は01:30、通知は01:43:23にlocal sendmailが受付。受信箱到達確認とは区別し、今回追加メールは送信しない。

| 最終成功統計 | 値 |
|---|---:|
| 検証採用job | 10087 / 10225 |
| 採用PicoDst | 30261（138jobの414本を除外） |
| 元producer実測input coverage | 812,298,245 events |
| accepted cache全decode | 802,090,792 frames |
| QAで実読取したPico Event/Track | 4,407,486 events |
| He4–Λ信号SE | **1,022,830 pairs** |
| He3–Λ信号SE | 5,001,632 pairs |

He3AsHe4は同じHe3 populationであり独立に追加加算しない。今回はこのjob campaignの全成功統計を使用し、別枠local4標本を合算しない。2021 catalog全体の解析完了とは主張しない。除外内訳はQA rc2=66、rc3=19、source/local size不一致53、計138。成功＋失敗集合と元dispatch全列/対象coverage、通知gate、最終ROOT/PDF/PNG SHAを独立35検査で照合した。

- [最終確認図：fixed r=1.4m、proton/π、全域と原点zoom（1頁PDF）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/track_merging_fullstats.pdf)／[PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/track_merging_fullstats.png)
- [旧EventFieldと複数修正定義の比較（4頁PDF）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/track_merging_fullstats_comparison.pdf)
- [最終全96histのROOT](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/aggregate/qa_10087_20261111T162940Z.root)
- [作成・入力不変・最終数値](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/completed.json)／[独立検証](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/independent_campaign_validation.json)

最終ROOT SHA=`f0cb59e9ea4978e3e9b151fe086b37498d0036d2fc43df3cafc7c43a1f78e89b`。前回9272source/932826pairから90004pair（9.65%）増加。固定1.4m＋globalηの原点両軸±.0125はproton27/π73、±.025は138/283、±.05は863/1131。±.0125の原点square/隣接square ringの面積補正密度比は.730/1.043、±.025では.571/1.001。全統計でも鋭い原点過剰は見えず、πは周囲と同程度、protonは中心密度が低めというrawSEの所見は変わらない。merging効率やcutへの結論にはMEが必要。ClosestSpatialEtaはpairごとの半径選択という条件付けがあるためpeakの根拠に使わない。

全域40bins/zoom native80bins、smoothingなし、1頁図はpanel独立color scale。raw信号SE・全保存k*/centrality、background subtractionなし、close-pair cut新規適用なし。入力ROOT不変、PDF1＋4頁・PNG・目視を確認。本番FillMergingのdebugとSE/ME診断の整備は上の未完了Todoを維持し、独立QAによる今回の統計拡張・図作成完了と区別する。


## 原点dip・β/曲率差の考察 — 07:49 JST / 2026-11-12

ユーザーの「αとprotonのβ/運動量/曲率差でtrackが重ならないのか」に対し、既存最終ROOTだけでdip形状とfrozen QA実装を再確認。job・hist・source/config変更は無し。[数値・説明用計算・一次資料](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/dip_interpretation_20261112.json)。

曲率半径はrhoT=pT/(0.299792458*|Z|*|B|)で、β単独から重なりは判定できない。同じrigidityならmassの違いでβは違っても曲率は同じ。説明用pz=0,B=.5Tの例では、proton pT1GeV/β.729とα pT2GeV/β.473はともにrho6.67m。逆に同じβ・方向ならmass/Zの違いでαのrhoはprotonの約1.986倍。この数値は実データの運動量分布ではなく、NIST massを使った反例計算。別曲率でも初期方向次第で1.4m一点に交差でき、他の半径では離れる。曲率差が長い並走を減らす可能性と、原点領域を禁止するという主張を区別する。低k*のα–Λは親の速度が近く、protonはΛ boostを受けるため、全sampleの平均β差を近接pairのβ差と同一視しない。

mergingは2trackの分離不成功によるpair inefficiencyなので、dipになり得る。splitting等は人工的な近接pair過剰を作り得るが、merging QAでpeakが必須ではない。[STAR二粒子解析](https://hepnp.ihep.ac.cn/article/doi/10.1088/1674-1137/abdf3f)と[ALICE一次資料 §3.3](https://arxiv.org/pdf/1507.06842)を参照。曲率と半径は[PDG2025 detector review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-particle-detectors-accel.pdf)の式で確認。他実験のcut値を採用したわけではない。

fixed1.4m/GlobalEtaのorigin halfwidth.00625/.0125/.025/.05のproton count7/27/138/863、中心/隣接square-ring面積補正密度比1.05/.73/.57/.82。πは22/73/283/1131、1.29/1.04/1.00/.97。最小4binだけが空になる鋭いholeではなく、数×10^-2程度の広い低密度領域。SpatialEtaでもproton8/24/120/859、比1.50/.75/.49/.82でdipが残り、GlobalMinusAsinにも同方向の形がある。k*<.2でも原点近傍pairは残るが少数。ratioはraw countの記述でME比や効率ではない。

実装では角度による原点mask/close-pair cutなし。fixed GlobalEtaのΔφ*はその半径上の位置φ差、縦軸はgMom方向η差なので同じ3D位置とは限らない。SpatialEta版は位置z差を反映する。StPicoTrack helix.originはPVへのDCA参照点でΛの物理的崩壊vertexではなく、実hitが存在する区間は未検査。7半径からpair毎に最短距離を選ぶClosestはadmissionも変わり、fixedの効率基準に使えない。

残Todo（未実施、今回新jobなし）：

- [ ] 同じ選択・幾何学でSE/MEのΔφ*–Δηを比較。cent/Vz/Bとsecondaryのorigin/PV整合を取り、rawSEの基準密度と同一event固有の欠損を分ける。SE/MEだけで物理相関を除去した絶対pair効率とは主張しない。
- [ ] 原点近傍のpairについて、物理β差と符号付き逆rigidity差Zα/pTα−Zp/pTpを保存・比較する。現exportには個別pairのβ/rigidity履歴なし。
- [ ] Λ崩壊以後の実際に共通するTPC padrow/hit区間で距離・近接継続区間/割合を調べ、1.4m一点の交差と長い並走を区別。近接効率の確定にはtruth/embedding等が必要。
- [ ] k*/rigidity依存・PID/track qualityを分け、cut採用や高k*CFの原因への結論は検証後に行う。


## 追加確認の再処理範囲 — 08:07 JST / 2026-11-12

userの再投入要否に対し、1)2D SE/ME、2)inverse-rigidity差、3)helix予測近接区間は、保存candidate＋元Pico Event/TrackのQA専用再読で実施可能。KF/PID・Λ候補選択の全面再実行は不要。ただし現ヒストにME/個別pairのrigidity/半径履歴は無く、candidateにも娘global helix幾何が無いため、この情報は現在のROOTだけから復元できない。全統計を調べるなら対象分のPicoを再読する追加QA jobが必要。

Pico StPicoTrackにはgMom/DCA origin/charge/helix、nHitsFit等、topologyMap/iTpcTopologyMapがある。topologyMapのbitはpadrow占有で、同じrowにhitを持つことと同じclusterを共有することは異なる。StPicoArraysの全branchにTPC各hitのXYZ/cluster ID/track-hit対応は無く、実hit共有には上流データ、失われたtrackのpair効率にはtruth/embedding等が必要。前節の「padrow距離」はまずhelixからの予測距離として区別する。

次の処理は1〜3を同じQA jobへまとめ、必要なselected trackのhelix/gMom/origin/charge/row mapとevent PV/B・provenanceを軽量geometry cacheへ保存する設計が有用。その後のME/曲率/半径依存はcacheから試行でき、Pico再転送の繰り返しを避ける。小標本は動作/速度の検証に使い、dipの判定には十分なpair統計が必要。今回source変更や追加job投入は行っていない。


## ユーザー指示による追加track-merging確認の保留 — 09:37 JST / 2026-11-12

ユーザーの指示により、前節で提案した追加2D SE/ME、符号付き逆rigidity差、helix/geometry cacheのためのPico再読は保留とする。これらは未実施のままで、新しいQA jobは投入しない。再開はユーザーの指示を待つ。既存10087成功sourceの最終ROOT・図・検証・完了通知、および原点dipの考察結果は保持する。本番FillMergingの修正・merging効率の確定・cut採用が完了したという意味ではない。

今回の作業は、CFを割る前の1次元k*のSEとMEを0.4–0.5 GeV/cで高さ合わせして重ね描きする確認へ移る。[CF高k* Todoの作図進行](cf_highk_systematics_20261111.md)を参照。この1次元SE/MEは既存1CEの物理ヒストから作図し、保留中の2D track-merging用MEとは区別する。


## SE/ME k*重ね描きの作図完了、追加merging QAの保留維持 — 09:50 JST / 2026-11-12

ユーザー依頼の1次元SE/ME重ね描きは、既存1CEヒストから0–20%・20–80%・0–80%の図を作成し、主検証と保存ROOT再読まで完了した。[主図0–20% PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_cent0_20.png)／[3頁PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_se_me_overlay_20261112_001/se_me_overlay.pdf)／[CF高k* Todoの検証記録](cf_highk_systematics_20261111.md)。0.4–0.5 GeV/cの積分を合わせたMEは約0.6 GeV/cから高k*へ向かってSEより上になる。raw分布形状差の確認であり、SE/MEどちらの歪みか、mergingが原因かを確定していない。

この作図はtrack-merging用2D ME・逆rigidity・helix/geometry cacheのためのPico再読を実施したものではない。これらの追加確認は引き続きユーザー指示で保留し、未完了Todoと既存track-merging全成功統計の成果物を保持する。新job・背景差引き・純度補正・本番source/config変更・mail送信なし。


## 解析コードdebug完了 — 01:49 JST / 2026-11-13

ユーザーの今回の明示指示でproduction FillMerging修正を実装し、batch-matched build・独立helix回帰・実Picoの修正hist/aliasを検証した。[詳細](../../../../analysisnote/20261113/lambda_he4_epd_ylab_merging_implementation.md)。旧physics SE/ME1675hist・約4000万cellsは不変。新close-pair cutなし。旧全成功統計の図を更新したものではなく、元から保留中の2D ME/逆rigidity/geometry cache・効率調査は引き続き未実施。本番再投入・mailなし。


## 2026-11-17：旧候補cache廃止後の傾向確認Todo

ユーザーが旧2AF候補ROOT13,027本の削除を明示指示した。作成済みの最終PDF/集約QA ROOT/検証・source記録を保持する。[削除と再確認のTodo](../17/old2AF_retirement_and_merging_recheck.md)。

- [ ] 新しい`anaFemtoLambda_4He_tree`からΔφ*–Δηを再作成し、`final_fullstats_20261112_004/track_merging_fullstats.pdf`と、proton/π⁻の全域・原点zoomの傾向が一致するか確認する。fixed r=1.4m/globalηなどの定義とevent/PID/統計差を照合して比較する。

これは今後実施する再確認であり、新しいQA job投入や比較完了を意味しない。

旧2AFの13,027本は2026-11-17に削除完了・残存0本を確認した。既存全統計のPDF/集約QA ROOT/検証資料は指紋不変で保持。上の新treeからの傾向比較は未実施のTodoとして維持する。
