# 未回収34件の原因解明（保留）

**状態：保留・未解決。** 2026-11-19のユーザー指示により、未回収34件の追加調査・再投入は当面行わず、今後の原因解明課題として残す。原因解決済み・損失が無視できるという扱いにはしない。第二段階は回収済み4,024 producerをまとめた20期間のROOTを対象に準備し、この調査をその前提条件にしない。

## 対象と判明済みの事実

- 解析：2021 P24iy、`anaFemtoLambda_4He_tree`。mainconf：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml`。
- SUMS ID：`9037A6F05FAD08FC8E9BD9286C8A38A5`。4,058 jobs中、正常回収4,024、未回収34。
- 内訳：unreadable PicoDst 32件、empty PicoDst 1件（subjob794）、chain Make途中失敗1件（subjob3413）。外側のSUMS exit0は解析成功を意味しない。
- 前回EP1/EP8と全worker入力groupが一致し、同じ34件で同分類の失敗が再現。32 unreadable＋1 emptyは同じ問題Pico LFNまで特定済み。3413の原因Picoは未特定。
- 未回収jobに属する入力は612/73,054 Pico（0.837736%）。うち310本はEntries実測あり（309 positive＋1 zero）、計8,594,396 Entries、残る302本は不明。既知入力欠損割合の下限0.464210%はaccepted event・Λ/He候補・pairの欠損率ではない。
- 明示的な原因Pico33本以外の579本は正常完走の証明がない。run別file欠損には最大36.69%等の偏りがある。
- run22162008の唯一のPicoは実測0 Entriesであり、正のevent統計をrun全体で失ったという意味ではない。zero-entryである理由は別途調べる。

## 再開時の調査項目

- [ ] 32 unreadableについて、元Pico/replicaの問題、転送・staging、アクセス権・一時的なIO障害を切り分ける。同一LFNで再現するだけで恒久的なデータ破損とは断定しない。
- [ ] subjob794のempty Picoについて、元データでも空なのか、取得・転送したファイルの異常なのか、既存ログと入力実体を照合する。
- [ ] subjob3413の失敗位置を入力listとevent位置へ対応付け、原因Pico・失敗branch/entryを特定する。同じSTAR/設定で最小入力による再現を検討する。requested502,094 / processed501,031 / accepted494,141 / status4で、未処理差1,063だけを損失と扱わない。
- [ ] 問題入力とその他入力を切り分け、579本を含む未確認入力の健全性を評価する。正常に読める入力だけを分離して回収できるか、再処理の要否を判断する。
- [ ] 回復を試みる段階で入力・成功cohortを明記し、accepted/candidate/pair統計およびrun偏りへの影響を測る。現在の20期間ROOTとの重複を避ける回収・統合方法を決める。
- [ ] 原因・再現条件・対処・回収範囲を記録してから、本課題の完了を判断する。保留中に自動で再投入しない。

## 再開に使う保存資料

- [日本語監査結果](../../../../analysisnote/20261118/tree_missing_inputs_20261118.md)：run別欠損・Entries下限・未解明点。
- [個人submodule内の証跡](../../../../verification/tree_run_period_hadd_20261118_001/README.md)：`failure_coverage.json`、`failed_worker_inputs.tsv`、`run_coverage.tsv`、`current_prior_comparison.tsv`。
- 親study原本：`rootfile/femto_lambda_tree_grouped_merge_20261118_001/failure_audit/`（監査scriptと上記4資料）。各失敗job・入力LFN・元ログの対応を保存。
- 初段状態記録：`rootfile/femto_lambda_tree_production_20261117_001/status_20261118_1018.json`と`tree/notify/status.json`。
- 実行時設定・native SUMS入力list・stdout/stderrは上記監査資料と保存runmetaから追跡する。原因調査に必要な失敗ログ・入力対応・runtime設定は、ログ整理時も保持する。
- [run時期別haddと後段の残作業](../18/tree_run_period_hadd.md)。未較正81runの保留とは別の課題として扱う。
