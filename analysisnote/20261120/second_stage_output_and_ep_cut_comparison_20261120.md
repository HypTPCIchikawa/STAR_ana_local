# 第二段階ROOTの保存内容とEP・pair cut比較

[3頁日本語PDF](../../../share/figure/20261120_lambda_second_stage/lambda_second_stage_guide_ja.pdf)。

既存depth100・EP8・cutなしの20出力を対象に、sourceと代表subjob0 ROOTのkey/branch metadataだけを確認した。event entries・histogram payloadは再走査していない。

## 実保存内容

20出力のlogical bytes合計は 11,491,533,504。代表1本のhistは1311個、監査treeは2種類、metadataは40個。

| histogram class | 個数 |
|---|---:|
| TH1F | 282 |
| TH1I | 1 |
| TH2D | 835 |
| TH2F | 169 |
| TH3D | 24 |

root直下204・true/576・mix/531。eventLedgerのcycleは18/19、selectedLambdaは1/2。論理的なtreeは2種類で、同じtreeを二重加算しない。

### eventLedger

91,040,515 entries。枝：`inputIndex`, `runId`, `eventId`, `status`, `cent9`, `selectedLambda`, `sourceInputIndex`, `candidateStoreStream`。

### selectedLambda

1,914,793 entries。枝：`inputIndex`, `runId`, `eventId`, `protonIndex`, `pionIndex`, `mass`, `px`, `py`, `pz`。

これらはeventと単体Λの監査台帳で、全pair treeではない。第一段階rich treeのfemtoLambdaReducedEvents/FemtoLambdaReducedMetadata、He4/娘track bank、全pairの相手・kstar・混合履歴を第二段階に複製していない。規格化済みCF keyは無く、raw SE/MEを別後処理で射影・規格化する。

## 比較条件

同じ20period、同じEPvalid母集団、depth100、bufferAll/bothDirections、元producer境界reset=falseを固定する。provider/較正/保存済8bin labelを保持し、後段poolとpair histogramだけをN1/8/16/32へ再binする。EP1は分割なし。

cut有無を比較する。ユーザー確定のveto対象はΛのproton娘–物理He4だけ。r=1.4m・実event BのSTAR helix、娘global momentum etaを使い、abs(dphi*)<0.05 rad AND abs(deta)<0.05を同時に満たすpairをSE/ME同条件で除外する。fullmass/signal/左右sidebandへ同じ受理判定を適用する。event全体・候補・pool入場を捨てず、πをveto条件に使わずQAを保持し、He3/He3AsHe4 controlへcutを掛けない。

通常の半径未到達pairは保持し、専用unreachable counterへ記録。必要raw bankの欠落・識別不一致はデータ不整合としてfail-loudで停止する。

新条件codeは編集中。実データ検証・新jobは未実行。cut前/後のsignal SE/ME 2D QAと12分類counterは追加予定で、既存1311 histには含まれない。

EP pool・mixing depth・pair vetoを変える場合、旧ME相手集合やpair間geometryを保存hist/audit ledgerだけから再構成できないため、第一段階rich treeから再mixが必要。規格化域や保存済軸のrebinだけなら後処理で対応できる。compact保存範囲を越えて候補選択を緩める場合は別問題として扱う。

## 証跡と範囲

- [metadata inventory](../../../share/figure/20261120_lambda_second_stage/second_stage_inventory_v2.json)
- [全key一覧](../../../share/figure/20261120_lambda_second_stage/second_stage_keys.csv)
- [全audit branch一覧](../../../share/figure/20261120_lambda_second_stage/audit_tree_branches.csv)
- [metadata集計](../../../share/figure/20261120_lambda_second_stage/metadata_summary.json)
- [比較設計](../../../share/figure/20261120_lambda_second_stage/comparison_plan.json)

共有production source/codeを編集せず、新規doc helperと資料だけを作成した。元ROOTの全event再解析・新規投入・メール送信は行わない。初段未回収34workerと未較正81runは保留を維持する。


[変更前source snapshot](../../../rootfile/femto_lambda_ep_merging_systematics_20261120_001/source_before_manifest.json)。
