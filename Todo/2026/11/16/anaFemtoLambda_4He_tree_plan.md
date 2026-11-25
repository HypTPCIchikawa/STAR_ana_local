# anaFemtoLambda_4He_tree 実装計画

状態：2026-11-17に実装と10,000-event検証を完了し、batch build・保存ROOTからのEP8/EP1/pairY5再解析とclosureがPASS。同日にユーザー指示で旧2AF ROOT13,027本を削除し、soft余裕2,095.44 GiBを確保した。代表sample・本番分割の評価は残す。本番job投入・commit/pushは未実施。[実装・検証結果](../17/anaFemtoLambda_4He_tree_implementation.md)。

2026-11-17更新：ユーザー回答により、第1段のROOTへ必要なQAヒストとtreeを同梱する方針に確定した。全ヒストをtreeだけから再計算することは要求しない。目的は、EPやpair rapidity等の変更時にPicoDst再読・再構成の本番jobを流し直さず、保存ROOTから再解析できることである。

同日追加決定：全統計1,842,342,504 eventと実測容量を照合し、ユーザーは省容量設定を標準にすることを承認した。`storePrePidControl: false`を標準とし、最終PIDを通らないcontrol候補だけを省く。初段108QA、selected-pair PID/rigidity QA、通常CF/EP/rapidity/mT・mergingに必要な候補/娘/eventは保持する。control全保存版は`fullControl`の別設定として残す。

## 1. 目的と処理の分担

PicoDstを一度読んで、Λ・³He・⁴HeとΛの娘p/π⁻、必要なevent情報をtreeへ保存し、診断に必要なヒストも同じROOTへ残す。そのROOTを後から何度でも読み直し、EP・centrality・rapidity・rigidity・mT・mixing条件を変更してSE/MEを作る。後段の再計算は必要だが、保存範囲内の条件変更でPicoDst再処理を繰り返す必要はなくす。

```text
第1段：PicoDst
  → 既存event選択・centrality・EPD Ψ1・KF Λ再構成・核PID
  → 同じROOTにevent単位のreduced tree＋必要な非pair QAヒスト
     （SE/MEを計算しない）

第2段：保存ROOTのリスト / TChain
  → 同じ共通pair計算 → SE・ME・pair診断・変更条件でのQAヒスト

第3段：ROOT後処理
  → SE/MEの比較・normalization・CF・sideband/purity処理・図
```

MEを第1段で作らない。第1段はtree専用出力ではなく、必要なヒストを残す再構成段階である。第2段ではPicoDstを開かず、KF再構成も再実行しない。CFの規格化・background/purity処理は既存方針どおりMakerの外に置く。

- 新しい入口の名前は`anaFemtoLambda_4He_tree`。runnerは`run_anaFemtoLambda_4He_tree.C`。
- 基準は11/13本番のΛ–⁴He＋³He control機能。既存のhistogram解析とcandidateStoreは互換性を維持する。
- 入力は2021 P24iyのみ。初期の比較対象は較正済み388run。未較正81runは引き続き保留し、既存係数の流用やrun上限の延長で有効扱いにしない。
- この変更で高k*のCFが1へ近づくとはまだ判断しない。新方式は原因検討のためにmixingと選択をPicoDst再読から切り離す。

## 2. 現行実装・oura方式との関係

| 対象 | 確認できたこと | 今回の採用方針 |
|---|---|---|
| 現行`StFemtoMaker`のcandidateStore | accepted eventとcanonical Float候補を保存し、同じpair関数でreplayできる。ただしwriterでもSE/MEを計算する | 型・長さ・参照・provenance検証、Float値の復元を参考にする。SE/MEを行わないtree＋QA producerを新設する |
| 同storeの保存項目 | 娘IDはあるが、娘momentum・raw dEdx・詳細TOF・helix・Λ fit情報が不足 | 既存storeを有効化するだけでは要求を満たさない。新しいschemaを明示する |
| 同storeのfile制約 | ROOT UUID変化でhadd/repackを拒否し、元producer境界を要求する | 新treeの永続的segment IDと、輸送用ROOTファイルのIDを分離する |
| `StCandidateFemtoMaker/CandidateStore` | 再構成とpairの2段方式。全track QAはproducer側ヒストを後段へ加算する | parts/persistenceの分離とQA持越しを参考にする。今回は第1段ROOTにQAを保持し、後段出力へ複製する場合は重複と母集団を管理する |
| ouraの`StFemtoPhiTreeMaker` | event tree・track treeをUIDで対応付け、後段が保存ROOTからSE→ME→pool登録を行う | 保存後の再解析、元track参照、保存範囲の明示、streaming読出しを参考にする |
| ouraの保存型・EP・pool | 主に整数/scalar branch、確認したEPはψ₂。候補ゼロのaccepted eventもpoolへ登録する | TVector3・今回のEPD Ψ₁・現行Λのpool入場条件を独自に保証する。整数packingとpool条件はそのまま移植しない |

ouraの確認対象は保存済みsourceであり、最新production設定との対応は未確定。`StPiMaker`という名前の実体は今回の調査範囲では見つからなかった。

## 3. 推奨するtreeの構造

### event内でtrackを共有する

基本案は、1 event = 1 entryの`events` TTreeとする。event情報と、可変長の候補・track bankを同じentryに置く。Λの娘が複数候補で共有されても、物理trackのraw情報は1回だけ保存する。

```text
events[1 entry]
  event情報・selection flags・元処理順
  track bank  ← Λ娘p/π⁻、選択済み³He/⁴He候補の和集合
  lambda候補  → protonTrackSlot / pionTrackSlot
  he3候補     → trackSlot
  he4候補     → trackSlot
  optional prePID control候補 ← fullControl設定のみ、必要なcanonical P4/PIDを保存
```

`lambdaP`・`he3P`・`he4P`とtrack bankのmomentumを`std::vector<TVector3>`として保存する案を第一候補とする。momentumはGeV/c、energyはGeV、massはGeV/c²、vertex/originはcm、BはkGとしてbranch仕様に明示する。ROOT5/SL24yで辞書生成・書出し・再読・圧縮量を確認してから最終branch形式を確定する。ouraのscalar branch実績を、その検証の代わりにはしない。型名だけTVector3と呼んで実体を別形式にすることはしない。

候補とtrackを別TTreeへ分ける場合は、絶対entry番号/offsetに依存しないUID joinが必要になる。初期案ではevent内参照を用い、この問題と全event→track mapのメモリ増加を避ける。

### 保存項目案

| 区分 | 保存する情報 | 目的・注意 |
|---|---|---|
| event識別 | dataset/production ID、runId、eventId、source file ID/entry、producer/segment UID、元input ordinal | 物理eventの重複検出と、元jobの境界・順序の再現を分けて行う |
| event選択 | 各段階の通過bit・棄却理由、再構成実施flag、元candidate数、pool入場情報 | event行の保存とpoolへの登録を区別する |
| centrality | cent9/cent16、refMult/fxtMult、補正multiplicity、centrality weight、nTracks/nBTOFMatch | 現行event QAとcentrality別の再解析 |
| vertex・磁場 | PVのTVector3、使用するvertex情報、B（kG） | Vz分類・DCA・helixの再現 |
| EPD Ψ₁ | 補正済みΨ₁、raw/recentered Qx/Qy、weight、hit数、valid/status、run/calibration index、元EP8 bin、元mixBin | 丸めによるbin変更を防ぐ。較正/定義ID・hashも記録 |
| track共通 | 元Pico配列indexとtrack ID、charge符号、role/species仮説bit、quality/selection flags | track ID単独で別eventのtrackを同一視しない |
| track運動量 | raw primary momentum、raw global momentumを別TVector3で保存 | 核の測定rigidityと物理momentum、primary/globalを混同しない |
| 飛跡 | origin（cm）、DCA、hits fit/max/dEdx、helix再構成に必要なstate | 必要に応じcurvature/dip/phase/helicityも保持。STAR版の実装と交点一致で十分性を確認 |
| TPC PID | accessorのdEdx（keV/cm）、nσ(e/π/K/p)、使用する核nσ/pull、較正定義、best-species判定 | raw ADCとは区別する。Pico nσ・核較正pull・KF PID pullを別名にする |
| TOF | hit/traits存在・match品質・値のvalidity、β、使用するtime/path length・matching量、m²等 | 現行selection/QAが参照する量を保存。TOFなしをβ=0の有効測定として扱わない |
| 核候補 | species/核のZ、trackSlot、physical momentumのTVector3、元pair用Float px/py/pz/E・charge・必要なkinematic cache、selected/prePID flags | trackのcharge符号と核のZを区別し、元candidate値は保持する。HeのZ=2補正は1回だけ。³He_as_⁴Heは同じtrackを参照しmass仮説だけ変更する |
| Λ候補 | momentumのTVector3、元pair用Float px/py/pz/E・cache、measured massとerror、娘slots、reco/selection flags | signalだけに限定せず、現行providerが返すfull mass cohortを保持しSB情報を落とさない |
| Λ fit | decay vertex、χ²/NDF・topology χ²/NDF、daughter distance、PV distance、decay/line lengthとerror/significance、cos pointing | 娘のmomentumだけではfit/topology QAを再構成できない |
| Λ娘 | track bank中のp/π⁻のmomentum・飛跡・PID/TOFに加え、候補ごとのDCA/path length・fit関連量 | global/primary momentumと、decay vertexで評価したfit momentumを区別する。後者を保存する場合は定義を明記 |
| event集計 | proton/pion仮説数・finderのΛ数等 | `hNProtonCand`等を無関係な粒子列なしで再現する |

初期版はlosslessな型を優先する。既存候補のcanonical Float値はそのまま保持し、TVector3からenergy/y/η等を再計算してbaselineと同一とみなさない。特にΛのpair energyは、measured invariant massではなく設定のpair massから作られている。

Heのpair momentumは測定primary momentumの2倍だが、merging用のglobal momentum/chargeへZ=2を二重適用してはいけない。m²についてもrigidity²から計算した量か物理momentum²からの量かをbranch定義に記載する。

EP invalidやTOF欠損には明示的なflagを持たせ、元のsentinel値をvalidな物理量として使用しない。EPD hit全体・PicoDst全track・raw detector hit・KF再fit用の全covarianceは基本保存対象にしない。raw Qの保存だけで、任意のEPD hit選択や未較正81runの較正が完了するわけではない。

## 4. 採用方針：第1段ROOTにtreeと必要なQAヒストを同梱

11/17のユーザー回答により、旧A/B比較のうち「必要な対象情報をtreeへ保存し、全track・cut前等のQAは第1段のヒストとして残す」方針を採用する。全既存ヒストをtreeだけから再計算するために全quality-trackを保存する案、QAのbin寄与を別TTreeに置き換える案は標準案にしない。

第1段でall-track/PID前QAを集計するためにPicoのtrackを参照しても、それだけを理由に無関係な粒子を永続treeへ残さない。診断用ヒストと、条件変更に使うevent/candidate/娘・必要controlのtreeに役割を分ける。

### ヒストの保存先とtreeの役割

| 対象 | 第1段ROOTに保存するヒスト | 後段tree解析で行うこと |
|---|---|---|
| event選択・centrality・vertex | cutflow、必要なcut前/後のVz・multiplicity・cent・weight QA | 保存eventの条件別集計。初段で棄却したevent全体の候補再構成は行わない |
| EPD Ψ₁ | status/validity、Ψ₁分布、cent別Ψ₁・元EP区分等の必要QA | EP bin数・pool条件を変更し、各planeのSE/MEを再計算する |
| Λ・娘・核の非pair QA | Λ mass/topology/acceptance、必要なdEdx・nσ・TOF/m²等の1D/2D QA | 保存候補と娘の情報から、変更条件で必要なQAを別に作る |
| 全quality-track・PID前QA | 必要なdEdx/TOF・核仮説QAを初段でfillして保持 | 無関係trackを再保存して全体QAを再生成することは要求しない |
| SE/ME k*・mass/mT/rapidity等のpair QA・CF入力 | 初段では作成しない | EP・centrality・rapidity・rigidity/mTの新条件をSE/ME両方へ適用する |
| merging Δη–Δφ* | pair依存なので後段へ置く | 保存したraw geometry・Bから交点/Δφ*/Δηを再計算する |
| prePID pair QA | 非pairのPID前ヒストは標準でも残す | 標準compactでは後段のprePID pair QAは対象外。fullControlのみTPC/best選択前bankから作る。selected-pair PID/rigidity QAは両設定で維持 |

必要なhistはYAMLの明示リストで選ぶ。event/cutflow、cent/vertex/EP、Λ mass/topology/娘acceptance、核・娘のPID/TOFを優先し、既存ヒストを無条件に全bookingしない。対象ごとにROOT key・型・軸/bin数・fill箇所・選択段階・母集団・weight・保存段階・必要branch・容量見積りを一覧化する。

`hNuclearPidPairJoint`等の巨大な多次元pair/sparseや、使わないchannel・重複した細分化ヒストは初段へ戻さない。THnSparseであることだけでは小容量とみなさない。1D/2Dでもbin数が大きければ容量/メモリを測定し、削減する際はQAの用途と比較対象を明記する。weighted histに必要なsumw2は保持する。

### 初段の固定QAと、後段の変更条件QAを区別する

第1段QAは保存時のevent/track/PID選択で集計した母集団を表す。後段でpair rapidityやEP binを変更しても、その初段ヒストが変更後の条件に自動的に対応するわけではない。例えばall-track dEdxヒストを「特定pair rapidityのpairに使われたtrackのdEdx」として表示しない。条件別の診断が必要な場合は、保存範囲内の対象候補・娘を使って後段で別にfillする。

初段QAは元ROOTに一度だけ保持し、後段の条件ごとの出力へ無条件に複製しない。レポートは初段QAと後段pair出力を参照して作れる。QAをまとめたROOTも必要になった場合は、producer/segment/QA集計範囲のUIDを記録し、同じQAの重複加算を防ぐ。treeの一部だけを再解析した出力へ、初段全体のQAを部分sampleのQAとして加えない。

QAの完了状態・母数・selection/weight・binning/hash・集計したsegment/entry範囲をmetadataへ残す。初段の固定QAと後段条件別QAはdirectory/名前/metadataで区別し、同じkeyを両段でfillして足し込まない。

省容量標準では、最終PIDを通らないcontrol候補を保存せず、変更条件でのprePID pair QAを再計算の対象外とする。selected-pair nσ・rigidity QAは選択済み候補だけから再計算できるため保持する。必要な場合だけfullControl設定でcontrol候補のcanonical P4/PID/選択flagsを保存し、条件別prePID pair QAを可能にする。初段のsingle-track集計ヒストを、変更条件のprePID pair QAの代わりとして扱わない。

この能力は必須設定`storePrePidControl`、ROOT metadata、保存時fingerprintへ明示する。compact/fullControlを誤った設定で読み込んだ場合は停止する。control-only trackへ不要なgeometryを複製せず、共有track・He仮説の重複も避ける。

cut前/棄却eventのQAは初段ヒストで保持する。これらを再現する目的だけで、全棄却eventの詳細や無関係trackをtreeへ追加しない。初段の通常event選択を通って再構成したeventは、EP valid/invalid・候補ゼロを含めて保存する。`baselinePoolEligible=false`を理由にevent行や候補を落とさず、元のpool入場条件・順序・counter等を保持する。これにより、保存時のEP適用やbin選択を後段条件へ固定しない。棄却eventには未実行の再構成結果や有効でないcent/EPを捏造しない。保存時より緩いevent cutで未保存候補を復元することはできない。

## 5. 後段MEの再現と系統検討

### baseline再現モード

新方式の検証時は、元job/segment境界・event順・候補順・候補仮説・pair前の選択を維持し、SE→ME→pool登録の順を固定する。元のVz/cent9/EP bin、FIFO、buffer深さ、maxMixEvents、両方向、pair weight/capを再現する。

現在のΛ解析は、物理candidate bankのどれかに候補があればpoolへ入れる。d/tだけのeventもpoolを占有し得る。d/tのmomentumは省略できても、`hasAnyLegacyCandidate`と元のselection/EP条件を含む`baselinePoolEligible`、元候補数を保存する。対象候補が空のeventも、このflagに従って空の対象payloadとして同じ位置でpoolへ入れる。

Λのみのevent、Heのみのeventも落とさない。最終SE pairが存在するeventだけへ絞る変更は行わない。prePID bankの有無を物理candidateのpool入場条件へ勝手に置き換えない。

乱数でpair capを適用する設定では、RNGの種・実装・draw順を検証契約に含める。seedだけで再現できるかは他処理のRNG消費を調べ、必要なら独立したmixing RNG/stateまたは検証用draw記録で比較する。未検証のrandom-cap結果までbitwise一致を主張しない。

### 条件を変える再mixモード

baselineが一致した後、同じtreeからEP1/EP8、EP bin数、Vz幅、buffer深さ、event gap、mixing方向等をYAMLで変えられるようにする。保存時のreconstruction/PID envelopeのprovenanceと、後段のmixing provenanceを分離する。後段でmixing条件を変えるたびにproducerのmixing hash一致を要求する設計にはしない。

初段の保存候補をEP8の特定planeや現在のpair rapidity範囲で間引かない。保存した連続Ψ₁と候補P4から後段で分類・pair選択する。EP status/validityは保持し、baselineの入場再現情報と、新しい後段mixing条件による入場判定を区別する。

- EP8の各planeのSE/ME k*を後段で個別出力する。MEの両eventに同じEP/Vz/cent分類を要求する。
- EP1でplane別表示を追加する場合は、MEのどちらのeventのΨ₁で分類するか（例：Λ側event）を明示し、EP8の同一pool bin表示と区別する。
- ⁴He YLabの-0.5分割、pair rapidityの現行10分割、rigidity/mTは、物理P4からSE/ME両方へ同じ定義を適用する。境界の等号処理・範囲はYAML/既存実装に合わせて検証する。
- rapidityでpairを選ぶことと、rapidityをeventのpool軸に追加することは別設定にする。片方の粒子の保存時にpartner未定のpair rapidity cutを掛けない。
- 保存時より緩いPID・topologyや、保存していない粒子の新channelは復元できない。保存envelope外の要求は明示的に停止する。

## 6. merge・読出し・provenance

- 初期運用は元treeファイルの明示リストをTChain/streamingで読む。不要なhaddの二重保存を避ける。
- ROOT fileのUUIDとは別の`producerUID/segmentUID/sourceOrdinal`を保存し、tree分割境界と元jobのpool reset境界を混同しない。
- `metadata` TTreeにschema、branch定義/単位、dataset、input manifest、処理順、元segment、完了状態、source revision/差分hash、実効YAML、centrality/dEdx/EP較正ID・checksumを保存する。segmentごとの予定/処理件数も照合する。
- 第1段QAのkey一覧・選択定義・weight・binning/hash・完了状態・集計範囲/UIDも記録し、split/repack/haddや後段出力の合成で同じQAを二重加算しない。
- 読出し時はschema・型・配列長・候補/娘参照・UID重複・ordinal順序・完了状態を確認する。不整合を黙ってskipしない。
- hadd/repackは新schemaの設計上の検討対象とするが、metadataの結合、segment重複、順序、分割を含むclosureに合格するまで対応済みとはしない。
- 元境界でpoolをresetするbaselineと、複数producer streamを連続してmixする解析モードを明示的に分ける。後者の統計/形状変化は系統として比較する。
- 初期版のhadd/repack検証は、完了したproducer segmentを丸ごと保つ順序付き結合を対象とする。同じsegmentを複数fileへ分割する形式は対応済みとしない。source input manifestは元producerにも保持し、hadd後も元ファイル/manifestを保管する。

## 7. コードへの反映順序

新しい専用pairing Makerは作らず、現在のΛで使用している共有`StFemtoMaker`を使う方針とする。新しい名前は解析entryの名前であり、Λ専用の別pair実装を意味しない。

| 手順 | 主な変更候補 | 内容 |
|---|---|---|
| 1 | 保存schema/track view/IO helper（`include/`・`StMaker/common/`） | 値だけのpayload、単位・型・参照とreader/writer契約。KF/Picoのobject pointerを永続化しない |
| 2 | `FemtoLambdaLegacy`・Λ providerのbridge | 再構成結果とraw track bankを取得できる共通処理に整理。PID・fitを二重実装しない |
| 3 | `StFemtoMaker`のΛ producer経路 | tree＋必要QAのproducer mode。明示リストの非pair QAとtreeを同じROOTへ保存する。pair関数とpool登録は呼ばず、巨大pair histもbookingしない |
| 4 | 同Makerのtree replay経路 | 新schemaを値のviewへ復元し、既存SE/ME関数へ渡す。Pico依存の娘acceptance/mergingはtrack viewから呼べる共通helperへ整理 |
| 5 | 新runner/analysis/run script | `run_anaFemtoLambda_4He_tree.C`、`anaFemtoLambda_4He_tree.C`、対応script。StChain Init/Make/Finishを使用し、`libStCommon`を先にload/link |
| 6 | producer/pair別mainconfと参照YAML | mode・input/output・schema・保存species/roles・storePrePidControl・保存envelope・QA明示リスト/軸・圧縮設定・後段mixingを明示。未知key/欠落/混在入力で停止 |
| 7 | 出力回収と説明 | joblistのtree＋QA同梱ROOT回収、QA/処理完了metadata、component README・REFERENCE・job/run手順を更新 |

MakerにKF libraryへの新しい直接依存を作らない。providerの注入と既存再構成partsを維持する。buildルールは既存の自動検出を使い、TVector3/STL辞書に必要な規則は既存方式を確認して追加する。

共通コード/再利用手順は共有mainへ、個別設定・study用ROOTマクロ・ログ・成果物・この計画/検証記録は`HypTPCIchikawa/STAR_ana_local`へ分ける。今回は計画ファイルの保存までで、新しいcommit/pushは実施しない。

## 8. 検証の受け入れ条件

| 検証 | 合格条件 |
|---|---|
| ROOT I/O | SL24y/ROOT5環境でTVector3と全branchがround-trip。欠損TOF・候補ゼロ・共有娘・複数候補を含め値と参照が一致。必要QAが同じ初段ROOTへ保存される |
| 既存physicsの維持 | 共通処理の整理前後で候補数・候補値・event選択・元one-pass histが一致 |
| 2段処理closure | 同じPico sampleで、one-passの対象QAは初段保存QAと、pair histはtree→pair出力と比較する。対応表の対象ごとにbin content、sumw2/error、entries、under/overflowと必要なhist statsが一致 |
| QA母集団・重複防止 | 初段固定QAと後段条件別QAを区別する。複数条件出力・split/repack/haddで初段QAが重複加算されず、部分tree処理に全体QAを誤って対応付けない |
| mixing | 元segmentごとのpool入場数・占有・候補数・pair counts・両方向・boundary resetが一致。d/t-only・Λ-only・He-only・全候補空を含む |
| merging | 元helixからの最初の前方cylinder交点、到達可否、Δφ*・Δηとhistが一致。radius m→cm、TPC z制限、charge/Bを確認 |
| EP/rapidity | EP8各binと全体の和、EP境界、invalid EP、He YLab=-0.5、pair Y10の境界、SE/MEの同一cut適用が一致 |
| input/metadata | 重複event/file、異schema・異dataset・不足較正、broken娘参照、途中終了、順序不整合を検出して停止 |
| TChain/merge | 元ファイル列と、許可したsplit/repack/hadd列で同じbaselineを再現。未検証の形式は拒否 |
| quota | 下記の実測projectionが容量/inode上限と運用余裕に収まり、全統計投入時の最大同時使用量が説明できる |

全統計jobの前にbatchと同じSTAR toolchainでbuildし、小さいround-trip、代表Pico sample、複数file/segmentの順に検証する。実装前の現在は、これらの新tree検証は未実施である。

## 9. quotaと計算時間の評価

新schemaの圧縮後サイズは未測定。過去の10,000-event pilotはhistogram ROOTが4,057,009 bytesという測定であり、今回のtree容量や全統計容量へ外挿しない。

実装後、central/peripheral・候補多重度・TOF有無を含む代表sampleで、compact標準とfullControlそれぞれのtree＋選択したQAヒストを測る。QAを追加/削減する比較では、そのリストとbinningを記録する。現時点では同一10,000-event sampleでfullControl約23.8 MB、compact約3.10 MBを確認している。全統計の容量保証には代表性とproducer単位のQA占有を別途評価する。

- 入力event数、保存event数、Λ/He候補数、重複除去後track数、prePID/QA track数。
- branch別の圧縮bytes、圧縮後bytes/event、bytes/候補、bytes/control track、読み書き時間。
- 初段QAのkey別圧縮サイズ、tree/QA/metadataそれぞれの容量、booking/fill時のpeak memory。bin数や選択QA数を変更した場合も実測する。
- ROOT logical bytesとGPFS allocated bytes、file/inode数、実際に適用されるquota。

tree部分のprojectionは`N入力event × 実測tree bytes/event`をcentrality等の代表層ごとに合成する。ヒスト部分は選択した軸/bin数・占有率・圧縮を含む実測値で、producerファイル数も考慮して別に足す。小sampleのhist容量をevent数へ単純比例させない。file overheadと高多重度sampleの余裕を含め、候補-only、prePID付き、異なるQA構成の見積りを混ぜない。

最大同時使用量は「既存ROOT + 新producer ROOT（tree＋初段QA）+ 残すsubjob + 後段hist + 一時書込み/コピー + haddする場合の重複出力」で評価する。同じproducer ROOTをtree分とQA分で二重計上せず、実fileのallocated bytesで検算する。イベント数と保存件数は2021 P24iy/較正済み入力の実測から求め、catalog file数だけを全統計event数と呼ばない。

pair rowやME相手の組合せは保存しない。再現できる派生量・共有娘・He仮説を不要に重複させず、圧縮/AutoFlushをYAMLで指定する。lossy整数packingや精度削減はclosureと容量測定後の別検討とする。

quotaの現在値とprojectionは実装pilot後・本番投入直前に取り直す。この計画段階では、新方式がquota内に収まるとも、全jobの所要時間が何時間とも断定しない。既存データの削除はこの計画に含めない。

## 10. 実装前に確定する事項とTodo

- [x] oura方式・現行store・必要hist・保存項目を調査。
- [x] MEを保存ROOTの後段へ移す方針、schema案、再現/容量検証を文書化。
- [x] 11/17のユーザー回答により、最小限の必要tree＋初段ROOT内の必要QAヒストを採用。全track/PID前QAはヒストで残し、これだけのために無関係trackをtreeへ保存しない。
- [x] 容量測定後の追加回答により、prePID control省略のcompactを標準に確定。selected-pair PID/rigidity QAと初段108QAは維持し、fullControlも別設定で保持。
- [x] 実装時にSTARの当該`StPicoTrack::helix`を確認し、必要stateとTVector3辞書を確定する。
- [x] 初段に残す必要QAの明示リスト・軸/bin数・容量見積りと、hist↔保存段階/branch対応表・保存envelopeを作成する。
- [x] 共通処理整理・tree＋QA producer・readerと後段pairを実装する。
- [x] ROOT5 round-trip、physics/merging/mixing closure、代表sample容量評価を通す。
- [ ] 単一10k試行の外挿433.35 GiBに対し、旧2AF削除でsoft余裕2,095.44 GiBは確保済み。複数run/centralityの代表性と実分割・QA/一時領域の評価後に本番jobを準備する。
- [ ] EP有用性の評価後に検討する未較正81runの較正は、引き続き別Todoとして保留。

### 根拠の場所

- 現行event→reco→store→pairs：`StMaker/StFemtoMaker/StFemtoMaker.cxx:350`、`:408`、`:432`、`:3152`。
- writer/replayの範囲と制約：`StMaker/StFemtoMaker/README.md:191`、`include/FemtoCandidateEventStore.h:25`、`:40`。
- canonical Floatとtrack参照：`include/FemtoCandidate.h:80`、`:130`。
- Λ fit情報/event集計：`include/FemtoLambdaProvider.h:15`、`:39`。
- 核rigidity・PID/control・acceptance/merging：`StMaker/common/FemtoLambdaLegacy.cxx:769`、`:804`、`:868`、`:976`、`:1128`。
- 既存2段方式とQA持越し：`StMaker/StCandidateFemtoMaker/README.md:143`。
- oura schema：`/gpfs01/star/pwg/oura/star-analyzer/include/FemtoPhiTreeSchema.h:317`、`:547`。
- oura producer：`/gpfs01/star/pwg/oura/star-analyzer/StMaker/StFemtoPhiTreeMaker/StFemtoPhiTreeMaker.cxx:447`、`:722`、`:1241`。
- oura後段：`/gpfs01/star/pwg/oura/star-analyzer/job/merge/phiqa_full_20260919/runtime/tools/anaFemtoPhiTreeDownstreamV3.C:1041`、`:1379`、`:1665`。
- 過去のhist容量pilot：`rootfile/femto_lambda_pair_rapidity_validation_20261113_001/quota_pilot_storage.json`。

[作成日の解析ノート](../../../../analysisnote/20261116/summary20261116.md) / [11/17の更新記録](../../../../analysisnote/20261117/summary20261117.md) / [Git保存先訂正の記録](git_remote_correction.md)
