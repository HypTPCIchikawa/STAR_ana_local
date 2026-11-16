# EP・YLab・merging実装の検証

2026-11-13 JST完了。ユーザー依頼の解析コード反映とローカル検証。本番SUMS再投入なし。

- [実装仕様・新3 mainconf・較正coverageと残Todo](../../analysisnote/20261113/lambda_he4_epd_ylab_merging_implementation.md)
- `build_command.json`, `build.log`：SL24y / gcc4.8.5 / ROOT5.34.38 batch-matched build。
- `epd_unit_retry.log`：EP20105チェックPASS（初回はtest loaderの不要な-lを修正）。
- `pair_unit.log`：YLab/helix462910チェックPASS。
- `cache_unit.log`, `cache_result.json`：schema1/2の138チェックPASS。
- `noEp_closure.json`, `noEp_closure_final.log`：1675physics hist・39963446cells完全一致、Y分割19392cells閉包、He4SE9/ME251不変。修正merging proton7/pion9、全alias一致。
- `larger_epd_mixing.json`：10000events、accepted9761、EPvalid9374/rejected387、SE853hist・20001030cells完全一致、He4SE13/ME459→403。He3SE58/ME2497→2111。
- `epd_mixing.json`：最初の1992events pilot。He4SE0だったため別runの1万eventsを追加した。
- `epd_run_coverage.json`, `epd_*runs.tsv`：全投入snapshot80982filesの90.21017%、388/469runs対応。未対応81runは追加較正が必要。
- `*_provenance.json`, `*_command.json`, `pico_*.log`, `larger_pico_*.log`：入力URI/SHA、catalog event記載とactual Pico entries、明示commandとStChain完了記録。
- `pico_noEp_qa.pdf`：16頁、全頁render・目視確認は`noEp_pdf_review.json`。
- `pico_epd8_qa.pdf`：17頁、EP頁renderは`epd_pdf_check.json`、`epd_qa_page17.png`。
- `implementation_manifest.json`：最終source/config/test/library SHA。
- `changes_from_preexisting_worktree.json`：編集前404filesのbackupとの比較。元からdirtyだった内容を保存。

`checkImplementation.C`の途中の失敗は検証macroのhistキーと保存directoryおよびimplicit Poisson varianceの参照を修正した。最終のclosureだけを完了判定に用いる。過去のfailureログは診断経緯として保存している。実装後に本番physics histが変化したという結果ではない。

小統計QAの空norm域ではCF未定義を明示する。今回の試験は高k* CF問題の解消やmerging効率を確定するものではない。元の追加2D ME・逆rigidity・geometry cacheの調査は保留のまま。
