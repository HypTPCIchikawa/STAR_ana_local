# run period groupingの最小証跡

20群の元producer件数・source run範囲・ROOT bytes・成功log event数、split runのgroup別Pico本数、10/20/30比較を保存する。ROOT本体・PicoDst・full snapshotは含めない。

元の詳細は `rootfile/femto_lambda_tree_grouped_merge_20261118_001/grouping_analysis/computed_final_run_cohorts/`。最終snapshotのSHAは `source_snapshot_final.sha256`、コピー前後のSHAは `manifest.json` を参照する。1 MiBを超えるコピーにはgzipを使う。

[日本語の選択理由と全20群表](../../../analysisnote/20261118/tree_run_period_grouping_20261118.md)。本証跡はgrouping計画を記録し、hadd全20群の完了やstage2実行を意味しない。
