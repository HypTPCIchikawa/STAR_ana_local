# 旧2AF候補ROOT 13,027本の廃止記録

2026-11-17のユーザー明示指示により削除した。対象は`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/`直下の完全prefix `auau3p85_anaFemtoLambda_4He_KFParticle_2AF731F2C75544890B4E858354DDE43D_`＋数値subjob＋`.root`に限る。全て通常file、本人owner、nlink1。

- `targets.jsonl.gz`：削除前13,027件のname/subjob/inode/device/size/blocks/mtime/ctime。実データbackupではない。
- `prepared.json`：明示指示、manifest SHA256、保持10ファイルの指紋、その他直下22件の記録。
- `independent_precheck.json` / `independent_postcheck.json`：削除前後の独立read-only監査。
- `deleted_jobs.tsv`：manifest順序と削除済みsubjob、時刻。全13,027行。
- `result.json`：削除数、残存0、保持物指紋・その他直下不変を確認。
- `quota_before.txt` / `quota_after.txt`：GPFS user quotaの原文。
- `completion_summary.json`：結果とquotaの集約。

削除allocated合計2,072,761,335,808 bytes（約2.073 TB）。quota使用量の実測減少は2,072,428,412,928 bytes、削除後soft余裕は2,249.96 GB（2,095.44 GiB）。quota差にはsmall audit file追加やfilesystem反映時刻も含む。

[削除結果と未実施のΔφ*–Δη傾向比較Todo](../../Todo/2026/11/17/old2AF_retirement_and_merging_recheck.md)。比較元PDF/派生QA ROOT・1CE/EP1/EP8・newtree検証は保持。旧2AFそのものの全内容hash/backupは作らず、ユーザー指示に基づき廃止した。コード変更・新job投入・commit/pushなし。
