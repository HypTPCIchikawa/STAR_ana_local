# 新tree本番投入の証跡

[投入記録](../../Todo/2026/11/17/tree_production_20261117.md)。2021 P24iy、較正388run/73,054files、compact初段tree＋108QA。SUMS ID `9037A6F05FAD08FC8E9BD9286C8A38A5`、実分割4,058jobs。

ここにはbuild/smoke closure・入力固定・source/config/loadedlib指紋・一回限定launch設定・時点付きprogressを保存する。`operations/`は使用launcher/読取監査scriptのsnapshot。ROOT、13MB入力list、bulk stdout/stderrは含めない。

親`rootfile/femto_lambda_tree_production_20261117_001/`が稼働中の主記録。`tree/launch_state.json`、`submission_progress.json`、`post_submission_audit.json`を参照する。生成・投入中は全投入完了/通知監視稼働とは扱わず、native launcherが全入力coverageとrunmetaを検証後にcompletion-only監視を起動する。別tmux observerが保存資料5種・SHA・4,058contiguous入力list・73,054件coverageと通知process設定を自動監査する。

自動hadd/再投入は行わない。元producerとsource/input manifestsを保持。最終ROOTの全schema/coverage評価と後段SE/ME、Δφ*–Δη傾向再確認はjob完了後の作業。
