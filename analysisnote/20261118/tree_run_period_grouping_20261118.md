# Λ–⁴He tree：近接run periodによる20群化

## 対象と実行状態

2021 P24iyの較正済み388runを対象にした初段producer campaign `9037A6F05FAD08FC8E9BD9286C8A38A5` の成功ROOT **4,024本**を、近接するrun時期ごとの20群へまとめる。元producerのevent順・segment・ROOTを保持する。今回の処理はhaddまでで、後段pair/MEは起動しない。

本番tmuxは `he4_tree_period_hadd_20261118_001`、起動記録は **2026-11-18 18:09:57 JST**。18:12:55 JSTの記録では `group_00` のsource catalog 236/236を通過して実mergeへ進み、完了群は0/20。部分出力の生成と全群の検証済み完了を区別する。以後の状態は [merge_status.json](../../../rootfile/femto_lambda_tree_grouped_merge_20261118_001/merge_status.json) と各群の検証結果を参照する。

- 初段mainconf：[main_auau3p85_anaFemtoLambda_4He_tree.yaml](../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml)
- 投入時config：[config snapshot](../../../job/run/configlog/config_auau3p85_anaFemtoLambda_4He_tree_treeWrite20261117_9037A6F05FAD08FC8E9BD9286C8A38A5.txt)
- 作業証跡：`rootfile/femto_lambda_tree_grouped_merge_20261118_001/`
- 出力先：`rootfile/auau3p85_anaFemtoLambda_4He_tree/run_periods20261118_001/`
- 小さい集計証跡：[verification/grouping](../../verification/tree_run_period_hadd_20261118_001/grouping/)

## 固定した成功統計

最後のsubjob 2869は `requested=550829 completed=550829 reconstructed=543693 status=0` を確認した。旧4,023本snapshotを上書きせず、2869を追加した [source_snapshot_final.json](../../../rootfile/femto_lambda_tree_grouped_merge_20261118_001/source_snapshot_final.json) と [source_files_final.list](../../../rootfile/femto_lambda_tree_grouped_merge_20261118_001/source_files_final.list) を新規保存した。

全4,024本について、通常file・owner・bytes・allocated bytes・device/inode・mtimeを照合し、成功stdoutの最終markerとMaker完了counterを再読した。requestedとcompletedの一致、status 0、reconstructedとの一致を確認し、ROOTのstatを再確認した。全件のlog監査は [final_completion_log_audit.json](../../../rootfile/femto_lambda_tree_grouped_merge_20261118_001/grouping_analysis/final_completion_log_audit.json) に残す。

| 項目 | 値 |
|---|---:|
| 成功producer ROOT | 4,024、全て一度ずつ採用 |
| 成功source PicoDst | 72,442、LFN重複0 |
| source list中のrun | 387 |
| 成功logのinput events合計 | 1,842,806,823 |
| 成功logのreconstructed/accepted合計 | 1,819,447,887 |
| 元ROOT logical bytes | 471,489,115,031 |
| 元ROOT allocated bytes | 472,021,729,280 |

これらのevent数はproducer成功logの累計である。hadd入力ROOTのtree/segment/108QAの検査はdriver側の検証結果に従う。起動前の実2source smokeではtree・segment・108QAが一致し、差binは0だった。

最終snapshot SHA256：

```text
9502f9b8c09af154cc6232d654c4026e0c836bab93f83f4b9662ec80c90fd2f0
```

旧4,023本snapshotのSHA256 `c1fd9caa2da2edcb6ea54856ddf75ea8b733057e0f5c6e27f71aa770ee694c0c` も一致した。未回収34workerを入力へ追加していない。

元の較正388runのうち成功sourceにrun全体が現れないのは **22162008**。このrunの元対象は `st_physics_22162008_raw_0000002.picoDst.root` 1本で、失敗logは `ERROR: empty PicoDst`。他runでも34失敗workerに属する計612 Picoの欠損が残るため、「387run」は欠損がこの1本だけという意味ではない。失敗内訳は unreadable PicoDst 32、empty PicoDst 1、chain Make failure 1。[欠損監査](../../../rootfile/femto_lambda_tree_grouped_merge_20261118_001/failure_audit/failure_coverage.json)を保持し、catalog全体の完全統計や欠損の物理的無視可能性は主張しない。未較正81runは保留のまま。

## run-min cohortと20群の選択

各producerのnative SUMS `.list` のfirst tokenから `st_physics_<run>_...` を読み、run集合・最小/最大runを求めた。ROOT payloadからrunIdを新たにdecodeした値ではない。元producerの内部順序は維持し、producer間だけを `(runMin, runMax, subjob)` 順に並べる。

4,024 producer中363本は複数runを含む。run区間の重なりを全て不可分にすると自然連結成分は18個だが、最大成分が **154,230,172,266 bytes** となり、20等分の目標約23.574 GBに比べて大きい。このため、同じ最小runを持つproducer全体を不可分にする **381個のrun-min cohort** を単位に、連続するcohortを20群へ分割した。個々のproducerやevent列は分割しない。

分割は容量の最大値を抑え、隣接境界の容量差を調整する。input/accepted eventsのバランスも別に評価した。20群のlogical bytesのCVは **1.60%**、input eventsは1.86%、accepted eventsは1.83%。各群185〜240 producerなので、全入力を同時openする上限は最大群240に8を加えた248とする。大きい連結成分の保持と、producer途中でのevent分割を行わず、近接runのまとまりと群容量を両立させる選択である。

| 案 | 1群のlogical GB（decimal） | 容量CV | split run | minority Pico | 全source Picoに対するminority |
|---:|---:|---:|---:|---:|---:|
| 10群 | 46.504〜47.772 | 0.94% | 8 | 74 | 0.1022% |
| **20群・採用** | **22.805〜24.236** | **1.60%** | **18** | **175** | **0.2416%** |
| 30群 | 14.458〜16.311 | 3.45% | 28 | 266 | 0.3672% |

10群は境界が少ない一方で1出力が約47 GBになる。30群は小さい出力になるが境界とrunの分散が増える。20群は約24 GB/群、統計量も近く、runの大部分を同じperiodへ保つ中間案として採用した。この選択からME統計やCF形状の改善を結論しない。

split runは、複数runを含むbridge producerを丸ごと保持した結果、同じrunのsource Picoが2群に存在するrunを数える。今回の18runはいずれも2群に限られる。runごとのPico本数が最大の群をmajorityとし、それ以外の群のPico本数をminorityと定義する。合計175本は全72,442本の0.2416%、split runに属する3,652本の4.7919%。これら175本もbridge producerとともに保持される。[run別・群別Pico本数CSV](../../verification/tree_run_period_hadd_20261118_001/grouping/groups_20_split_runs.csv)で追跡できる。

## 採用した全20群

「起点run範囲」は各群内producerの最小runの範囲。実sourceの最小runはその左端と等しく、bridgeが含む実source最大runだけを別欄へ示す。19境界中18境界に隣接runへのbridgeが残り、group15/16境界は分離している。ROOT bytesは元producerのlogical bytes合計で、hadd後サイズの保証値ではない。

| 群 | producer | 起点run範囲 | 実source最大run | ROOT bytes | input events | accepted events |
|---|---:|---|---:|---:|---:|---:|
| 00 | 236 | 22158032–22160027 | 22160028 | 23,438,571,158 | 92,891,006 | 91,040,515 |
| 01 | 214 | 22160028–22161006 | 22161010 | 22,805,353,302 | 88,485,719 | 87,421,760 |
| 02 | 228 | 22161010–22161030 | 22161031 | 22,963,912,798 | 88,641,747 | 87,596,220 |
| 03 | 240 | 22161031–22163010 | 22163011 | 23,813,610,691 | 91,863,814 | 90,818,334 |
| 04 | 195 | 22163011–22163029 | 22163030 | 23,600,912,694 | 92,186,007 | 91,046,915 |
| 05 | 198 | 22163030–22164008 | 22164009 | 23,531,594,298 | 91,639,189 | 90,505,160 |
| 06 | 193 | 22164009–22164030 | 22164031 | 24,235,987,508 | 94,838,244 | 93,673,595 |
| 07 | 194 | 22164031–22165014 | 22165015 | 23,114,038,352 | 90,023,092 | 88,960,957 |
| 08 | 195 | 22165015–22165034 | 22165035 | 23,541,495,762 | 91,853,411 | 90,741,888 |
| 09 | 194 | 22165035–22166029 | 22166030 | 23,508,757,919 | 91,912,136 | 90,789,339 |
| 10 | 198 | 22166030–22167010 | 22167011 | 23,606,846,433 | 92,370,673 | 91,257,378 |
| 11 | 191 | 22167011–22169002 | 22169003 | 23,486,094,782 | 92,017,129 | 90,867,761 |
| 12 | 191 | 22169003–22170007 | 22170008 | 23,509,508,516 | 91,994,052 | 90,866,756 |
| 13 | 197 | 22170008–22171030 | 22171031 | 23,827,330,004 | 93,361,703 | 92,194,024 |
| 14 | 192 | 22171031–22173002 | 22173003 | 24,007,737,517 | 94,142,319 | 92,970,453 |
| 15 | 201 | 22173003–22173023 | 22173023 | 24,142,231,004 | 94,660,400 | 93,455,508 |
| 16 | 198 | 22173024–22174046 | 22174047 | 23,093,046,797 | 90,691,018 | 89,528,409 |
| 17 | 190 | 22174047–22175008 | 22175010 | 23,651,981,749 | 92,800,132 | 91,589,405 |
| 18 | 194 | 22175010–22177005 | 22177006 | 24,137,540,191 | 94,502,852 | 93,329,021 |
| 19 | 185 | 22177006–22179022 | 22179022 | 23,472,563,556 | 91,932,180 | 90,794,489 |
| 合計 | **4,024** | | | **471,489,115,031** | **1,842,806,823** | **1,819,447,887** |

全4,024 subjobの一意性、2869の採用、20個のordered listとassignmentの一致を確認した。全割当と各input listは [groups_20/assignment.json](../../../rootfile/femto_lambda_tree_grouped_merge_20261118_001/grouping_analysis/computed_final_run_cohorts/groups_20/assignment.json)、10/20/30比較は [alternatives_summary.json](../../verification/tree_run_period_hadd_20261118_001/grouping/alternatives_summary.json) を参照する。

## 後段MEに必要な設定と未実行範囲

haddはproducer segmentのmetadataを保持する。既定のsegment境界pool resetを保ったままでは、20ファイルへまとめてもproducerを越えるMEは増えない。期間ごとにpoolを継続する後段は、専用の明示設定で次を揃え、**1 periodを1 job**として読む必要がある。

```yaml
resetPoolAtSegmentBoundary: false
requireOriginalMixing: false
```

この条件は元producer単位のmixingと異なる解析条件になる。元のstrict replay設定を上書きせず、period単位の設定・input manifest・統計・SE/ME比較を残す。event/cent/Vz/EPのmixing条件と保存候補の範囲は別に確認する。

現在は `stage2Automatic: false`、`stage2Started: false`。後段の本番投入・CF作成・MEの物理評価は未実行であり、自動起動しない。元producer ROOTを削除・移動しない。
