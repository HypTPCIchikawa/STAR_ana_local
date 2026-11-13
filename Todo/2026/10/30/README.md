# 2026-10-30

- [更新ROOTでのCF再計算・結果と検証](lambda_he4_cf_update_20261030.md)（完了）

- [当日のTodo・作業進捗](STAR_Todo20261030.txt)
- [Λ–4He CF：過去資料と更新ROOTでの再作成準備](lambda_he4_cf_preparation_20261030.md)
- [次回の入力・設定・実行前条件](next_run_plan.yaml)
- [当日の作業ログ](../../../../analysisnote/20261030/summary20261030.md)

## 過去の図付きレポート

- [総合レポート](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/summary_4He_KFParticle_3p85.pdf)（10頁）／[README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/README.md)
- [Λ質量窓の低側／高側比較](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves/mass_halves_comparison.pdf)（13頁）／[README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves/README.md)
- [mT別CF：0–10%](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10/cf_mt_cent0_10_constant_vs_method5.pdf)（8頁）／[README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10/README.md)
- [mT別CF：0–20%](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_20/cf_mt_cent0_20_constant_vs_method5.pdf)（8頁）／[README](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_20/README.md)

## 検証記録

- [更新ROOT検査](input_histogram_check.json)：175必要ヒストは全件存在・軸一致・非空。通常1,786ヒスト、event数一致。
- [実入力リストの照合](dataset_input_check.json)：81,082ファイルはすべて2021 P24iy、全26,994 subjobリストとの対応一致。
- [投入時設定との照合](submitted_config_check.json)：12設定の内容一致（2件はコメントのみ差）。
- [PDF一覧・頁数・hash](assets_manifest.json)／[13マクロの互換性照合](macro_compatibility.json)／[後処理settingsの保存コピー](settings/)

更新ROOTで13マクロ、総合10頁、質量窓比較13頁を完了。規格化fit不成立はCF全体をNAとし、hash・入力不変・欠測伝播を検証済み。

## 今回の完成PDF

- [更新ROOTの総合CFレポート（10頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/summary_4He_KFParticle_3p85.pdf)
- [Λ質量窓比較レポート（13頁）](../../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/mass_halves/mass_halves_comparison.pdf)
- [完成成果物・検証一覧](cf_outputs_manifest_20261030.json)
