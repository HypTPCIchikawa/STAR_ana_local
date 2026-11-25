# tree producer未回収34 jobの入力欠損監査

対象は2021 P24iy、3.85 GeV Au+Au固定標的の `anaFemtoLambda_4He_tree`。SUMS IDは `9037A6F05FAD08FC8E9BD9286C8A38A5`、mainconfは `config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml`。入力73,054 Picoを4,058 workerへ割り当て、正常回収4,024 ROOT、未回収34 ROOTを確認した。

監査は保存済みstdout/stderr、SUMS入力.list、runmeta、既存入力group署名と最終source snapshotを用いる。PicoDstやproduction ROOTの新規読取り・再解析は行っていない。本記録は未回収入力の監査であり、20 groupへのhaddの完了・検証結果とは独立する。

## 前回EP投入との一致

未回収34 jobの内訳はunreadable PicoDst 32、empty PicoDst 1（794）、chain Make失敗1（3413）。前回EP1 `FD44A5EB1B0AA22C546F169E86482587`、EP8 `DA0DEB4347BD834F1096F8D8798F2133` に対し、全4,058 workerのordered LFN group署名が同indexで一致した。今回失敗34の全入力URL・ordered LFN group・エラー分類も両campaignと一致する。

32 unreadableと1 emptyは同じ33 Pico LFNを特定できた。3413は同じworker入力groupとMake失敗が再現しているが、原因となったPico LFNは特定できない。unreadableという記録だけでは元Picoの恒久的破損と転送・アクセス失敗を区別できない。

今回失敗34は旧EPの共通除外35から3496を除いた集合で、新たな失敗workerはない。3496は旧EP1では正常ROOTがあったがEP8との共通selection外、旧EP8ではunreadableだった。今回は正常回収されたため、最終成功cohortは旧共通4,023に3496を加えた4,024となる。同job数だった途中の4,023-source snapshotは3496追加・2869未回収という異なるcohortだったため、最終snapshotで再監査した。

## 全体の欠損量と分母

| 項目 | 実測・集計値 |
|---|---:|
| 凍結入力Pico | 73,054 |
| 未回収workerに割り当てられたPico | 612（34×18） |
| 未回収Picoのfile数割合 | 0.837736% |
| 原因をLFNまで特定したPico | 33（32 unreadable＋1 empty） |
| 上記33以外の巻き添え入力 | 579 |
| preflightのEntries測定あり | 310（positive 309＋zero 1） |
| Entries未測定・不明 | 302 |
| 未回収入力の既知Entries合計 | 8,594,396 |
| 最終成功4,024 sourceのinput events | 1,842,806,823 |
| 最終成功4,024 sourceのaccepted events | 1,819,447,887 |

既知の未回収入力Entriesだけを用いた入力欠損割合の下限は、

`8,594,396 / (1,842,806,823 + 8,594,396) = 0.464210%`

となる。成功入力に対する既知未回収Entriesの比は0.466375%。302 PicoのEntriesが不明なので、0.464210%を最終欠損率の推定値や上限とは扱わない。

Entriesは各workerの `[input-file] entries=N` とPico LFNを対応付け、LFNで重複排除したTTree Entriesの実測である。全branch・全eventが正常に読み通せる保証ではなく、accepted events・候補数・SE/ME pair数の損失率でもない。失敗分accepted総数は未確定。SUMS .list/.datasetや既存入力metadataに利用可能な正のper-Pico event列はなく、catalogの未記録・placeholder 0を実測0 eventと混同していない。

579の巻き添え入力は「明示的な原因Picoとして記録されていない」という分類にすぎず、正常完走を証明したhealthy入力ではない。早期失敗により後続Picoが未チェックのままのworkerもある。309 positive Picoもpreflightでheader/Entriesを読めた範囲の確認である。今回、入力再分割・正常性の追加検査・リカバリー再投入は実施していない。

## runへの集中

612 Picoの未回収は28 runにまたがる。下表はfile数割合で大きいpositive/未確定入力のrunであり、event重み付き割合ではない。

| run | 未回収Pico / 全Pico | run内file数割合 | 既知未回収Entries下限 | Entries不明Pico |
|---|---:|---:|---:|---:|
| 22177012 | 51/139 | 36.69% | 559,207 | 26 |
| 22177013 | 39/113 | 34.51% | 800,592 | 11 |
| 22165006 | 72/239 | 30.13% | 812,834 | 45 |
| 22176001 | 53/203 | 26.11% | 712,295 | 27 |

既知未回収Entriesの最大runは22165006の812,834で、同runにはさらに45 Pico分のEntries不明が残る。

run 22162008は1/1 Picoが未回収でfile数上は100%だが、その唯一の `st_physics_22162008_raw_0000002.picoDst.root` はpreflightで `entries=0` と実測された。したがって、このrunを「正のevent統計が全損したrun」とは表現しない。同worker794にはrun 22162007の既知129,629 Entriesを持つ6 Picoと、後続run 22162009の未チェック11 Picoが同居し、worker全体の出力が回収されていない。

## worker3413の途中失敗

3413はrequested 502,094、processed 501,031、reconstructed/accepted 494,141、status=4。requestedはpreflightで得た全入力Entries、processedはその内数で、両者を加算していない。未処理の差は1,063 input eventsだが、正常なproducer ROOTとしてはworker全体が未回収であり、損失を1,063だけとは扱わない。原因Picoと全branchの健全性は未確定。

## 解釈と未完了事項

同じ33原因LFNと34 workerの再現は確認できたが、「統計損失が無視できる」とは断定できない。入力段階で少なくとも0.464210%に相当する既知Entriesが未回収で、302 Picoが未測定、run内file損失には36.69%等の集中がある。物理選択後のaccepted/candidate/pair損失、観測量のrun依存・mixingへの影響は未評価。リカバリーは未実施であり、可能性と実行済みを混同しない。

## 再現資料

root原本は `rootfile/femto_lambda_tree_grouped_merge_20261118_001/failure_audit/` に保持。解析scriptは同ディレクトリの `audit_missing_inputs.py`、入力はproducer controlの `status_20261118_1018.json` とgrouped studyの `source_snapshot_final.json`。再現コマンドはrepository rootから次を実行する。

```bash
python3 rootfile/femto_lambda_tree_grouped_merge_20261118_001/failure_audit/audit_missing_inputs.py \
  --source-snapshot rootfile/femto_lambda_tree_grouped_merge_20261118_001/source_snapshot_final.json
```

privateへ複製した4 artifactを以下に示す。サイズは非圧縮原本bytes、SHA256も非圧縮原本に対する値。1 MiB超の原本だけgzip複製とし、root原本は変更しない。

| artifact | 原本bytes | 原本SHA256 |
|---|---:|---|
| [failure_coverage.json](../../verification/tree_run_period_hadd_20261118_001/failure_coverage.json) | 940,634 | `4fb059be6d9f67150b764a173c9dd3c2a9969186bd55c237ee49340663bc9f8c` |
| [failed_worker_inputs.tsv](../../verification/tree_run_period_hadd_20261118_001/failed_worker_inputs.tsv) | 169,145 | `5299205622dfb3213a7f06ff94d6fe8fa4d11261f3ea2ab0477af755848efa9b` |
| [run_coverage.tsv](../../verification/tree_run_period_hadd_20261118_001/run_coverage.tsv) | 2,386 | `2accbad99f698d32b313e31ef3e179af5dcc2534ace8e63ef8036a06c3caa0c4` |
| [current_prior_comparison.tsv](../../verification/tree_run_period_hadd_20261118_001/current_prior_comparison.tsv) | 4,941 | `b691713d53f82bfd96d8bb18ad7e09e65874017edc369a5fa4028096eb06e8a2` |
