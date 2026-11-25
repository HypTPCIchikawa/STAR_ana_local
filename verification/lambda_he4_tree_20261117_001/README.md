# Λ–He4 reduced tree 実装の検証証跡

[日本語の実装・検証結果](../../Todo/2026/11/17/anaFemtoLambda_4He_tree_implementation.md)。初段QA＋treeと保存ROOTからの後段SE/MEを検証した。compactを標準とし、selected PID/rigidity QAは保持する。

`implementation_manifest.json`は使用source/configと外部較正/入力list・project librariesのSHA256を記録する。`source/`は小さいsource/configのsnapshot、`implementation_source.patch`は共有tracked sourceの差分。`study_macros/`は個別ROOT検証マクロ。生成ROOT・lib・bulk logsはここに含めず、親の`rootfile/femto_lambda_tree_implementation_20261117_001/`へ保持する。未較正81run・本番投入・commit/pushは対象外。

- compact QA108/pair1186と直接Pico、さらに変更前11/13 ROOTとの比較が全bin/error/entries/flow一致。
- EP1のSEは不変でMEは変化。pair rapidity10→5は再bin結果と一致。
- IO/完全segment hadd・metadata・capability・geometryとhist契約はPASS。
- 既存入口テスト2件は変更前snapshotでも同じ設定不整合を再現しており、全テストPASSとはしない。
- `storage_projection.json`は単一10k sampleからのtree外挿433.35 GiBとsoft余裕165.69 GiBを示す。全統計には保存先/空き容量と代表sampleの再評価が必要。

ROOT内source revisionはunrecordedであり、このmanifestを代替の完全なGit commitと解釈しない。元producerとordered input/source manifestsを保持する。hadd/repackは完全segmentの順序/metadataが保たれる範囲で検証した。
