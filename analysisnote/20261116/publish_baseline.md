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

---

## Session update — 16:00 JST / 2026-11-16：両remoteへの公開と境界確認

mainの内容公開commitは`00bf4bdd07d247a344bb126a6c2f84112706e099`、privateは`37d798a7428d89ddb6f0f86b66b7b2c5ecdf8cdc`。15:59 JSTに両GitHub remoteのrefs/heads/mainを読み、各ローカルHEADと一致、親のgitlinkもprivate commitと一致を確認した。公開差分は明示33path＋gitlinkだけで、個別設定・後処理macro・ROOT/ログ等の混入無し。

独立boundary監査でもmain33path、private575pathの分類がPASS。主要後処理source45blobが現在版と一致、退避9ファイルはactiveから消え、private drafts内だけに保存された。source/docs150ファイルのpublic indexは投入packageのSHAと一致し、生成.so/.d10ファイルは公開していない。private保存/復元/サイズ制約の10 testsもPASSした。

検証receipt・boundary監査・baseline SHA・build/test要約は[verification/publication_20261116_001](../../verification/publication_20261116_001/)へ保存する。この完了記録を含むprivate後続commitと親gitlink更新を追加するが、解析本体の内容は上記の公開commitから変えない。main作業treeに残る個別config/後処理overrideは私用snapshotで保存済みの意図した残差であり、一括stageしない。
