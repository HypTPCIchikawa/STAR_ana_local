# 旧2AF候補ROOTの削除とΔφ*–Δηの再確認

2026-11-17、ユーザーが旧2AF候補ROOT 13,027本の削除を明示指示した。比較用の既存PDF・集約QA ROOT・source/config/入力記録を保持し、後日新しいtreeからΔφ*–Δηを再作成して傾向を照合する。

## 削除対象と保持物

対象は`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/`直下の`auau3p85_anaFemtoLambda_4He_KFParticle_2AF731F2C75544890B4E858354DDE43D_<数値subjob>.root`だけ。13,027本、logical 2,071,420,278,951 bytes、allocated 2,072,761,335,808 bytes。旧候補cacheをバックアップせず廃止する指示であり、manifestは実データの代替ではない。

状態：**旧2AF 13,027本を削除済み、残存0本**（2026-11-17T14:24:17.602184+09:00）。削除allocated合計2,072,761,335,808 bytes。比較元PDF/QA ROOT等10件の指紋と、その他直下22件の不変を確認した。削除後soft余裕は**2,249.96 GB（2,095.44 GiB）**。対象manifest・journal・quota前後を保持する。

比較元は以下を保持する。

- [指定の最終確認図（PDF）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/track_merging_fullstats.pdf)
- [比較用4頁PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/final_fullstats_20261112_004/track_merging_fullstats_comparison.pdf)
- [最終QA集約ROOT](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_campaign_20261111_001/aggregate/qa_10087_20261111T162940Z.root)
- 上記studyの完了・独立検証記録、元source/configと入力manifest。1CE/EP1/EP8の結果、新tree検証も保持。

## Todo：新treeからの再作成と傾向比較（2026-11-20実施済み）

- [x] `anaFemtoLambda_4He_tree`で保存した情報を使い、Λ–⁴HeのΔφ*–Δηヒストを再作成する。Λ娘proton/π⁻をそれぞれ確認する。
- [x] 旧PDFと同じfixed radius r=1.4 m・global momentum ηの定義について、全域図と原点zoomを作る。必要な角度符号・charge/rigidity・飛跡の定義を旧QA実装と照合する。
- [x] 原点の鋭いpeakが見えない点、proton側で中心密度が低い傾向、π側が周囲と概ね同程度という傾向が再現するか、指定PDFと比較する。色スケール・軸・binningを揃え、中心/隣接領域の面積補正密度比も確認する。
- [x] event/run範囲、EP-valid選択、centrality・k*・Λ mass/PID・統計数の差を記録し、単純な総count差と形状差を区別する。新旧統計が異なる場合、全binの完全一致を合格条件にしない。
- [x] 比較PDFと数値・入力指紋を保存し、結果をこのTodoへ追記する。merging効率やcut採用の確定とは区別する。

今回の指示は再確認作業をTodoへ残すもので、この削除作業と同時に新jobやヒスト再生成を実行するものではない。旧2AFを元入力として使う過去の再実行手順は、削除後はそのまま利用できない。再確認は新tree経路で行う。

[削除manifest・保持物指紋・実行結果](../../../../verification/delete_old2AF_20261117_001/) / [容量調査](disk_capacity_audit.md) / [tree実装・検証](anaFemtoLambda_4He_tree_implementation.md) / [旧track-merging記録](../11/track_merging_audit_20261111.md)

## 2026-11-20：再確認結果

期間内ME深度100の全20出力のSE QAを集計し、proton側dip・π側ほぼflatの傾向を旧QAと照合した。中心±0.05/周囲の面積補正密度比はproton新0.816±0.021・旧0.821±0.032、π新0.967±0.023・旧0.972±0.033。新旧の母集団差を記録し、共通粗binの形状比較とnative-bin raw表示を分けた。exact cut件数/merging効率の確定ではない。

[比較PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/share/figure/20261120_lambda_ep_merging_systematics/track_merging_current_vs_old_20261120.pdf) / [EP・cut比較Todo](../20/ep_bins_and_proton_he4_merging_cut.md)

input stat/旧QA SHA・bare/EventField alias一致は検証済み。原本は`rootfile/femto_lambda_ep_merging_systematics_20261120_001/merging_plot_validation.json`。削除済み旧2AF ROOTは再使用していない。
