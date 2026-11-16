# 2026-10-31

- [追加依頼：Fig.13・14の質量分布比較と規格化域の点検（11/02に完了）](cf_normalization_followup_20261031.md)

追加確認では20–80%のMethod5の規格化を0.40–0.50 GeV/cへ移すと計算可能。ただし近隣区間での安定性は未検証。失敗fitの滑動窓には約203万カウントがあり、低統計とは説明できない。10/31は承認レビューtimeoutで停止したが、11/02に追加図と指定総合PDFの更新を完了した。

- [当日のTodo](STAR_Todo20261031.txt)
- [SE/ME純度fit点検とMethod3・5・6比較：結果・解釈・検証](method6_study_20261031.md)
- [明示設定](settings/method6_settings.yaml)
- [数値結果概要](method6_results_summary_20261031.json)

## 図・レポート

- [比較レポート（20頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/method6_cf_purity_report_20261031.pdf)
- [CF比較：2分割・3分割（2頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/cf_methods_3_5_6.pdf)
- [Method3 全fit（100頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/atlas_method3.pdf)
- [Method5 全fit（40頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_method6_20261031_001/pdf/atlas_method5.pdf)

## 記録

- [全145PDFの頁数・hashと完了記録](method6_outputs_20261031.json)
- [最終数値・来歴検証](method6_validation_20261031.json)
- [report metadata](method6_report_metadata_20261031.json)
- [ベクトル図の余白処理記録](method6_crop_manifest_20261031.json)

SE低統計と質量モデル不適合を確認し、Method6を指定どおり比較した。20–80%のMethod5/6はME規格化fit不成立によりNA。Method6を正解または最終採用とは認定していない。

- [個人snapshotの保存一覧](method6_private_snapshot_20261031.json)
- [本日の解析ログ](../../../../analysisnote/20261031/summary20261031.md)
