# 本番ベースラインの維持・main/private公開

## 15:43 JST / 2026-11-16：公開準備

ユーザーは「11/13本番投入時の解析機能を維持し、未実装のEP別k*下書きは退避してpush」を明示選択。mainは`star-tohoku/star-analyzer`、サブモジュールは`HypTPCIchikawa/star_local_submodule`、両branchはmain。fetch時点の両HEADはそれぞれorigin/mainと一致していた。

既存git HEADへresetすると本番機能まで失うため、reset/revertは行わない。EP混合・pairY10・He4単体YLab・merging診断を維持する。投入packageと現在のStMaker64/include40/analysis56、計160ファイルがSHA一致（analysisには生成.so/.d計10を含み、これらは公開しない）。既存config411本も投入時と完全一致、既存本番参照closureは新draft設定を参照しない。

未実装の設定3本・テスト2本・準備helper/fixture4本は[drafts/ep_kstar_20261116_001](../../drafts/ep_kstar_20261116_001/README.md)へ元パス/SHA付きで移動した。C++実装無し・未検証の下書きであり、activeなconfig/tests/rootfileから除外している。履歴や既存ROOTを削除しない。

mainには明示した共通実装・runner・必要header・job管理・対応検証等33ファイルを公開する。`common/macro/checkHistAnaFemtoLambda.C`、`tools/`、新しいROOT後処理macro・個別config・ROOT/log/build生成物はmainのstage対象に入れない。後処理macro6本の現内容は既存private snapshotにも保存済みであり、main作業treeの個別overrideは保持する。

privateは605pathをpermitted/size/sha検証してincremental snapshot更新（535copy、総記録22987、過去snapshot削除0）。新studyのsource/settings/小さな定量・来歴、33ページまとめと詳細PDF/図を保存する。121ページPDFだけは既存policyの個別例外機構で4MiB上限を明示し、一般上限1MiBは維持する。大きなROOT/so/圧縮jobbundle/ログ全量は保存対象外。`workspaces/`は従来どおりstageしない。

検証はbatch整合の`script/singularity_make.sh`（SL24y、sl73_x8664_gcc485、ROOT5.34/38）によるclean buildがreturn0。skills同期12source/12wrapper PASS。job管理testsはarchive20、notify16、merge_notify21、include-job選択9、合計66がPASS。monitor.lockのResourceWarningが一部mockテストに出たが、全test終了はOK。実データ処理・通知は起動していない。

保存した証拠は`rootfile/git_publish_baseline_20261116_001/`（`baseline_comparison.json`、`snapshot_before_draft_archive.json`、`draft_retirement.json`、`public_paths.json`、`publication_plan.json`、`private_capture_result.json`等）。private snapshotにsource/helper/監査JSONを保存。現時点はcommit/push前で、remote更新後に検証結果を追記する。

本番識別子はEP8 `DA0DEB4347BD834F1096F8D8798F2133` / EP1 `FD44A5EB1B0AA22C546F169E86482587`。mainconfは`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_epd{1,8}_pairY10_covered25.yaml`。今回の新規mainconf/jobid/configlog・本番job・Pico再読・メール無し。未較正81runの保留を維持し、EP別k*実装は別途の新指示まで再開しない。
