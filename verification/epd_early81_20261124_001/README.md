# 3.85 GeVの初期81 runとΛ–⁴He tree解析

- `retry_audit_20261125.json` に再投入54件の入力、終了コード、イベント件数、出力ヘッダーの確認結果を記録した。
- `retry_failure_log_excerpt.txt` に11番のROOT展開エラーと終了行を保存した。
- `campaign/` に較正・出力検証用マクロと初回復旧の記録を保存した。マクロに引数を渡す際は、元のcampaignのROOTファイルとリストを指定する。
- `snapshot/config/` にtreeの本番・比較設定と初期81 runのEP較正設定を保存した。初期81 runのrecenter・shift ROOTはSTAR上に保管している。`calibration_locations.tsv` の取得元とSHA256を使い、設定が参照する配置先に置く。再現時はanalysis infoの作業パスも確認する。
- 従来388 run用の較正ROOTはSTAR上の参照データを用いる。出典、取得元とチェックサムは `snapshot/config/event_plane/xialei_epd_2021P24iy/` に記載した。
- `snapshot/common/macro/` に出力ROOTのQAマクロを保存した。実行時はリポジトリの同じ相対パスに置く。
