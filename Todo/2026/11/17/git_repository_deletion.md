# 旧star_local_submodule repo本体の削除

2026-11-17、ユーザーから手動承認下でagentによる削除を明示依頼された。既存認証は見つからなかったため、公式CLIをchecksum確認して隔離した一時領域へ配置し、ユーザー本人のdevice OAuth認証完了後に実行した。

- [x] 削除対象を`HypTPCIchikawa/star_local_submodule`、repo ID `1429165992`に限定。
- [x] 認証アカウント`HypTPCIchikawa`と対象repoのadmin権限を確認。
- [x] 削除対象のrefsはmain `aaed3cd53a63aa1d420a0a5d141b89e61e46e6f6`のみ。旧個人HEADと元STAR_ana_local HEADの履歴が、保全先mainの祖先として残ることを再確認。
- [x] 09:30:51 JST、DELETE APIがHTTP 204で成功。
- [x] 09:32:46 JST、認証済GETがHTTP 404となり旧repo本体の消失を確認。
- [x] `HypTPCIchikawa/STAR_ana_local`（repo ID `1158941728`）が正常で、main `a128fcb31c4b6c075674db255763f25c1703707f`を維持。
- [x] `star-tohoku/star-analyzer`のmain `f35f4424f35a2bb8c06458721259ae21bf24c50d`と親gitlinkを維持。両remote HEADとも一致。
- [x] ローカル個人workspace・既存解析の手元差分を保全。
- [x] 隔離した一時認証ファイル（hosts.yml/config.yml）・一時CLIを消去。既存GitHub設定は変更しない。
- [x] 前日の削除Todoを完了へ更新し、[tree解析の計画](../16/anaFemtoLambda_4He_tree_plan.md)側へ戻れる状態を記録。

証跡は親解析ディレクトリの`rootfile/git_repository_deletion_20261117_001/`に保存した。`predelete_verification.json`、`delete_response.json`、`postdelete_verification.json`、`temporary_auth_cleanup.json`を参照。秘密値・アクセストークンは記録していない。要約は個人repoの`verification/git_repository_deletion_20261117_001.json`へ保管する。

この作業でGitのcommit/push、解析コード変更、PicoDst再読、job投入は行わない。tree解析は[既存の実装計画](../16/anaFemtoLambda_4He_tree_plan.md)の段階にあり、全track/PID前QAの保存範囲は未確定のままである。
