# EP共有補正・全run対応の可否調査

[結論と再現可能なsource根拠](../../analysisnote/20261113/ep_calibration_transfer_and_all_runs.md)。

- `review_metadata.json`: 調査完了・production47files不変・共有係数の適合性は未検証。
- `epd_eight_bins.png/.pdf/.root/.json`: 既存run22159051/10000event pilotの有効9374を全centで投影。左上は補正後Psi1とnominal45度境界、右上は実Float分類の0–7占有数、下段はcent依存。全統計の図や平坦性合格の証拠ではない。
- `DrawEpdEightBins.C`, `draw_command.json`, `draw.log`: 読み取り作図sourceと実行記録。
- `early_epd_presence/`: 前期run22123040先頭100eventのcompiled StMaker/StChainによるsignal存在確認。source/log/inputSHAとJSON。補正の流用適合性は確認していない。

Picoはgain補正済nMIPを保持するが、recenter/shift済EP角は保持しない。Centrality共通表のrun受付範囲拡張と、run別EP実測係数の追加は同じ操作ではない。共有係数を対象期間で独立検証できれば採用の選択肢になる。適合しなければEvent/Track/EpdHitによる軽量recenter→shift生成を行える。

未対応81runsはday121–125、現在の388runsはday158–179。隣接raw-Q QAのROOT509binsと現run表82rowsが一致しない候補は採用しない。新較正・production変更・SUMS投入は今回実施していない。
