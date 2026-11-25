# anaFemtoLambda_4He_tree 実装・検証結果

[確定計画](../16/anaFemtoLambda_4He_tree_plan.md)に従う実装と10,000-event検証を完了した。標準はユーザー承認済みの省容量設定。初段ROOTには必要QAとtreeを同梱し、SE/MEは保存ROOTを読む後段で計算する。本番job投入・commit/pushは行っていない。

## 保存範囲と再解析

新入口は`analysis/anaFemtoLambda_4He_tree.C`、runnerは`analysis/run_anaFemtoLambda_4He_tree.C`。既存StFemtoMakerの共通pair処理を使い、新しい専用pair Makerは作らない。

| 段階 | 内容 | 保存・再計算範囲 |
|---|---|---|
| 初段 | 2021 P24iy PicoDstのevent選択・EPD Ψ₁・KF Λ・核PID | event/cent/vertex/EP、Λと娘p/π⁻、選択済み³He/⁴Heの運動量・飛跡・PID/TOFをtreeへ保存。必要な非pair QA 108個も同じROOTへ保存。SE/MEは作らない |
| 後段 | 保存ROOTの再読 | 共通pair処理によるSE/ME、EP区分別k*、pair rapidity・⁴He rapidity・mT・rigidity・merging。PicoDst readerとΛ再構成は起動しない |
| 標準compact | `storePrePidControl: false` | 選択済みpairのPID/rigidity QA 540個は保持。非選択prePID controlを使うpair QA 540個だけを保存範囲外とする |
| fullControl | 同keyをtrueにした専用設定 | 非選択prePID controlも保存。compactとはcapability/fingerprintで区別し、未保存情報を要求する再解析を拒否する |

track bankは選択候補とΛ娘の和集合だけを保存し、娘を共有する候補間でraw trackを重複保存しない。候補のphysical momentumとraw primary/global momentumは実体が`TVector3`。pair計算用の元Float P4/cacheも保持する。Z=2はphysical momentumへ一度だけ適用し、merging用helixへ二重適用しない。

EP-invalid event、target候補が空のevent、d/tだけで元poolへ入ったeventも元処理順・入場flagとともに保存する。無関係なd/tの粒子列は保存せず、MEのpool占有を再現する情報を残す。後段で変更できるのは保存候補の範囲内の条件であり、候補再構成/PIDを緩めて未保存粒子を回復することはできない。

初段108QAは初段の固定母集団に対応する。後段の条件別QAへ自動複製・加算しない。CF規格化・purity/sideband処理はMaker外の後処理とする。高k*のCF形状が改善したとの物理結論は、この動作検証からは出さない。

## 実行入口

作業場所は親の`star-analyzer/`。全引数が必須で、mainconf/input/output/jobid/maxEventsを明示する。

```bash
bash script/singularity_make.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml --no-clean -j4
bash script/singularity_run_anaFemtoLambda_4He_tree.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml <PICO_LIST> <PRODUCER_ROOT> <JOBID> -1
bash script/singularity_run_anaFemtoLambda_4He_tree.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_replay.yaml <PRODUCER_ROOT_LIST> <REPLAY_ROOT> <JOBID> -1
```

`..._tree_replay_ep1.yaml`でEP1、`..._tree_replay_pairY5.yaml`でpair rapidity 5分割を作る。fullControlは`..._tree_fullControl.yaml`と`..._tree_replay_fullControl.yaml`。通常は元producerファイルを順序付きリスト/TChainで読む。producer境界でpoolをresetする標準が元jobとのclosure条件となる。境界を越えるmixingは明示的な別条件として扱う。

## 検証結果

実データは2021 P24iy・較正済みrun22159051の同一PicoDst先頭10,000 event。9,761 eventを再構成した。Λ174、³He2361、⁴He611、raw関連track3317を保存した。EP-invalid387、target空7282、d/t-only pool入場5526を含む。

| 確認 | 結果 |
|---|---|
| SL24y/ROOT5 batch整合build | PASS |
| TVector3/Float P4・rich Λ/TOFのIO round-trip、型/長さ/参照/metadata/重複/未完了拒否 | PASS |
| 完全producer segmentの順序付きhadd/repack | PASS |
| compact/full capability・選択fingerprint不一致の拒否 | PASS |
| native Pico helixと保存stateのgeometry、非空merging、娘acceptance等 | PASS、135,721 checks、20交点/10非交点 |
| compact初段QA108・後段pair1186と直接Pico解析の比較 | PASS、全bin/content/error/entries/flow差ゼロ |
| 11/13の変更前に作った同sample ROOTとの独立比較 | PASS、同じ108＋1186個の差ゼロ |
| EP区分別ヒストの総和32比較 | PASS、inclusiveと一致 |
| EP8→EP1 | 保存ROOTだけでPASS。SE599個は不変、MEは変更。He4 signalはSE13、ME403→459 |
| pair rapidity10→5 | 保存ROOTだけでPASS。8 TH3は10分割を2binずつ合成した結果と一致。他1194個は不変 |
| fullControlの初段108QA・pair1726、およびEP1/pairY5 | PASS |
| 既存hist契約 | 632比較PASS |
| 既存入口テスト | 9件中2件失敗。変更前HEAD＋元からの個人設定overrideを再構成したsnapshotでも同じ2件が失敗することを確認 |
| 既存ユーザーのtracked差分15ファイル | 作業前snapshotとbyte一致 |
| 入口script引数・YAML参照・writer/replay joblist生成 | PASS。XML生成のみで未投入 |

入口テストの既存失敗は`test_eight_mainconfs_reference_presets_and_real_macros`と`test_nuclear_pair_pid_qa_is_opt_in_for_canonical_production`。既存個人overrideの核PID設定参照の不整合であり、今回のために設定を変更したりテストを緩めたりしていない。最初のYAML省略型/ROOT5検証API/loaderの失敗ログは修正後PASSと区別して保持した。

## 容量・メモリと全統計への制約

| 同じ10,000入力event | compact | fullControl |
|---|---:|---:|
| 初段ROOT logical bytes | 3,100,258 | 23,832,870 |
| 初段ROOT GPFS allocated bytes | 3,145,728 | 検証JSON参照 |
| 保存TTreesの圧縮basket合計 bytes | 2,525,634 | 23,234,130 |
| 初段108QAの圧縮key合計 bytes | 520,198 | 520,198 |

compactのROOT全体は約1/7.69に減った。`time -v`で測ったcompact初段と後段のpeak RSSは301,536 KiB（約294.5 MiB）と402,336 KiB（約392.9 MiB）。wrapper/子processを含む10k試行であり、大きなproducer全体のメモリ保証ではない。

既存388run統計のeventLedgerとinputEvents parameterはともに**1,842,342,504**。単一runのtree bytes/input eventをこの数へ外挿すると、compactでも**約433.35 GiB**、fullControlは約4.28 TBになる。QAはevent数へ比例させず、実producer数・占有率とともに別途加算する。ROOT/file metadata、割当block、後段出力、一時コピーもこのtree外挿には含まない。

実装pilotの最終quota照会（旧2AF削除前）では、inDoubt差引後のsoft余裕は**165.69 GiB**、hard余裕は677.69 GiB。したがってcompactでも現soft余裕に収まらず、tree単体だけで約267.67 GiB不足する。さらにtreeをhaddで二重保持するとtreeだけでhard余裕も超える。**全統計の本番投入には、保存先/空き容量と実分割・代表sampleの容量計画が必要**。soft graceへの依存や、未承認のevent/候補削減・既存ROOT削除は行わない。

### 旧2AF削除後の容量更新（2026-11-17）

ユーザー指示で旧2AF候補ROOT13,027本を削除し、約2.073 TB分を廃止した。現在のsoft余裕は**2,249.96 GB（2,095.44 GiB）**。上の容量不足は削除前の評価であり、compact tree外挿433.35 GiBに対する基礎容量は確保できた。複数run/centralityの代表性・実分割とQA/一時コピーを含む本番容量計画は引き続き確認する。本番jobは未投入。[削除とmerging傾向再確認のTodo](old2AF_retirement_and_merging_recheck.md)。

## metadataと利用上の範囲

- treeは`femtoLambdaReducedEvents`、segment metadataは`FemtoLambdaReducedMetadata`。完全producer segmentの順序・entry数・物理event重複を検証する。segment途中の任意分割は未対応。
- hadd/repackは完全segmentなら検証済みだが、file-levelのsource TObjStringだけで全segmentの元source listを自己完結して保持するわけではない。元producerの入力manifestと順序を別途残す。
- `hasTof`はtraits存在、`tofValid`はβがfiniteかつ正のみ。全fieldの有効性やmatch品質通過を意味しない。β/time/match/localY/Z/hitPositionはnative値。現在native getterのない`tofPathLength`・`tofDeltaY/Z`は予約値-999で、値を復元したわけではない。
- embedded source revisionは`unrecorded`＋警告。汚れた作業ツリーはGit HEADだけでは再現できないため、source/config/dependency SHAと元inputをsidecarで保管する。
- 未較正81runは保留のまま。EP有用性を確認した後の独立較正・検証は別Todo。

## 証跡と残作業

親の`rootfile/femto_lambda_tree_implementation_20261117_001/`に実ROOT・command/log・closure JSON・quota・manifestを保存。個人repoの[検証証跡](../../../../verification/lambda_he4_tree_20261117_001/README.md)へ小さい検証結果とsource/config snapshotを保管する。generated ROOT・lib・bulk logは共有mainへcommitしない。

- [x] schema/Reader/Writer、共通Maker producer/replay、飛跡・娘再現、新入口・strict設定を実装。
- [x] compact/fullControl、初段QA108、後段selected PID/rigidity540、EP区分別k*を用意。
- [x] batch build、IO/geometry/metadata、実データclosure、EP1/pairY5、既存回帰との切分け、メモリ/容量測定を完了。
- [x] Todo・解析記録・実行説明を更新。既存差分とmain/private HEADを保持。
- [ ] 複数run/centralityの代表sampleで容量外挿を改善し、実分割・QA/一時領域を含む本番計画を確定する（旧2AF削除で基礎容量は確保済み）。
- [ ] 本番投入は別の依頼で行う。commit/pushも未実施。
- [ ] EP有用性の評価後、未較正81runのrecenter/shiftを同じ方法で準備・独立検証する（保留）。


### 2026-11-17追加Todo：旧merging図との傾向比較

- [ ] 新treeからΔφ*–Δηヒストを作り、[指定PDFを保持した再確認Todo](old2AF_retirement_and_merging_recheck.md)に従って、旧全成功統計のproton/π⁻・全域/原点zoomの傾向と比較する。旧2AF候補ROOTはユーザー指示で削除対象となり、比較元PDF/QA ROOTを保持する。
