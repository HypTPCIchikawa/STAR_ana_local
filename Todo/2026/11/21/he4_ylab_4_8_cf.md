# ⁴He単粒子YLabの4/8分割CF

ユーザー依頼：4He自身のラピディティを4分割ないし8分割し、高k*でのCFを調べる。pair rapidityとは別量として扱う。

- [x] 既存は−0.5の2分割のみで4/8復元不可。rich treeから第二段階だけ再実行する方針を確認。
- [x] producer20期間の固定QAを一度だけ集計。He4候補127,785,819、99.9104%が−1≤YLab<0。YLab<−1は114,452(.0896%)。全cent/EP-invalidも含む候補QAであり、CFのpair比率とは別。
- [x] EP8/cutなし、depth100、期間内継続、同20入力を基準として独立YAMLを作成。8bin[-1,0)・幅.125、4binは隣接2binのraw SE/ME・Sumw2加算後に規格化。
- [x] Legacyの単粒子He4YLab TH3をfull/signal/左右SB×SEMEへ追加。有限範囲外はUF/OFに残し、既存inclusive・2分割・pairY・EP・He3control・poolを維持する。
- [x] batch整合build・実Legacy境界/flow/ME方向/物理P4テスト・同400k保存eventで旧全hist不変/全Y投影closureを確認。
- [x] 同20期間の本番replayを一度投入し、quota/送付物/実入力/native資源を検証。
- [x] 4/8各Y領域のCF/scaledSEME/centrality比較、範囲外control、highk窓比を作図。完了後mail1通。

主比較はYLab[-1,0)、低Y tailを別枠で保持・報告し、全Y inclusiveの照合へUF/OFも含める。旧strict2分割はexact−.5を双方から除外、新8binはROOT半開区間で上側に含むことを境界検証で区別する。SE/MEは同じHe4候補自身のYを選び、event全体やΛ単粒子Yを切らず、mixing poolをYで分けない。

原本：`rootfile/lambda_he4_ylab_split_20261121_001/`。入口：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_he4YLab8_mix100.yaml`。元source snapshot24filesを保持。初段PicoDst再実行なし、未較正81run・未回収34producerの保留維持。

検証完了：batch build成功、実Legacy合成テスト93,775,373チェック（1,188 pair calls）成功。同じ保存event先頭400,000件で旧1,311ヒストの内容・誤差・軸・flow・Entries・Sumw2・統計量が完全一致、新8TH3の全Y投影も一致。最終pilot作図 `pilot_postprocess_003` は320曲線の検証・目視確認に成功。規格化窓の統計が足りない区分はNAとする。watch設定検証成功。本番は同じ20 rich ROOT（保存event合計1,819,447,887件）を一度だけreplayし、完了後の検証・4/8分割図作成・通知を行う。

証跡保全：`star_local/verification/he4_ylab_split_20261121_001/`。大きなROOT本体の複製なし。既存ソース24件中、変更は `FemtoLambdaLegacy.{h,cxx}` とmaker READMEの3件のみ、その他21件と既存設定・20入力の保全を確認。

本番確認済み：SUMS `FEDC5185F359546FB1B7B39F70DC8144`、20期間/20job。監視初回pollで20件running・terminal0。tmuxとwatch PID 3510711 の生存、startup pins、入力coverage、ZIP/runtime/設定/資源指紋を確認。送出後quota余裕1,128,065,662,976 bytes・grace none。全統計のCF結論は処理後に判定する。`handoff_verified.json` とprivate verification/handoffを保存。

---

## 2026-11-22 全統計の完了確認・結果

SUMS `FEDC5185F359546FB1B7B39F70DC8144` は20/20正常終了、保存event1,819,447,887件、EP-valid1,736,875,306件。20ROOT計11,494,216,306 bytesを回収。11/22 07:50:55 EST（NY）に成功mailを送信済み。ROOT再読のhist検証290,880箇所がPASS、全Y inclusiveは旧EP8/depth100/cutOff基準と一致。完了後の独立確認でも20ROOTのstat指紋、成果物12件のSHA256、7PDFの全98ページのデコードが成功、主要3ページを目視確認した。

条件は0–60%、EP8、pair cutなし、norm H=[0.4,0.5) GeV/c。下表は各窓でのSE/ME積分比をHで規格化したCFであり、bin CFの単純平均ではない。

| ⁴He YLab | CF[0.6,1.0) | CF[0.6,0.8) | CF[0.8,1.0) |
|---|---:|---:|---:|
| 全Y | 0.93023 | 0.95129 | 0.91216 |
| 範囲内合算[-1,0) | 0.93006 | 0.95112 | 0.91199 |
| [-1, -0.75) | 0.98943 | 0.98552 | 0.99319 |
| [-0.75, -0.5) | 0.97456 | 0.98332 | 0.96642 |
| [-0.5, -0.25) | 0.91947 | 0.94189 | 0.89975 |
| [-0.25, 0) | 0.89500 | 0.92706 | 0.87090 |

8分割のCF[0.6,1.0)（低Y側から）は0.98884, 0.98993, 0.98638, 0.96551, 0.93867, 0.90520, 0.89616, 0.88665。低Y側では1に近い領域がある一方、高Y側に高k不足・下向きの形状が残る。Yで選択しても全域の問題は解消していない。最も低い8binは規格化SEが1,178と少なく、形状の平坦性を強く主張しない。4分割は同じ8binの生SE/MEと分散を合算した結果で独立標本ではない。

今回のYはSE/ME双方のHe4候補に課すpair投影条件で、mixing poolはcent9×Vz×EP8・depth100のまま。背景差引き前のsignal-window CF。誤差表示はalpha固定のcount-onlyで、規格化窓の揺らぎや共有event/ME相関を含まないため有意度判定はしていない。0–60%のY<−1はSE4,105/2,221,644（0.1848%）、ME755,085/459,646,635（0.1643%）、Y>=0は0。これらは全Y inclusiveに保持される。

成果物の原本：`rootfile/lambda_he4_ylab_split_20261121_001/completion_watch_001/he4_y_results/`。`he4_y4_CF.pdf`, `he4_y8_CF.pdf`, 両`CF_highk.pdf`, 両`SE_ME_scaled.pdf` は各16ページ（0–60%、各10%区間、9つのcent9）、`he4_y_flow_controls.pdf` は2ページ。`cf_histograms.root` は320曲線の補助ROOT、`he4_y_cf_metrics.json` に数値・窓・raw countsを保持。主0–60%の12 nativeY曲線は規格化可能、ほかの空flow/低統計の条件にはNAを明示。

検証原本：`D/completion_review_20261122_001/{verification,numeric_summary}.json`。private保全：`star_local/verification/he4_ylab_split_20261121_001/completion_review/`。

残課題：高Y側に残るSE/ME形状差の原因を検討する。未回収34 producer・未較正81runの保留は継続。今回の確認で新job投入・source変更・再メールは行わない。
