# 個人Git保存先の訂正

ユーザー指定の個人保存先は[HypTPCIchikawa/STAR_ana_local](https://github.com/HypTPCIchikawa/STAR_ana_local)。`star_local_submodule`への公開は保存先の取り違えであり、今後のpush先として使用しない。

## 保全と訂正方針

- main訂正前：`79b36df8b8bf8db6a8194352e7598aa2fb6d8ac9`。
- 個人作業訂正前：`aaed3cd53a63aa1d420a0a5d141b89e61e46e6f6`。
- 正しい保存先の既存main：`50198217806dcc12725e3333feea9873e899801a`。
- 両既存履歴をmergeで保持する。正しい保存先の`Lambda_check/`・`Nuclear_id/`と、Todo・解析ノート・snapshot・報告・退避draftを保全する。
- ローカル配置`star_local/`とTodo/analysisnoteのsymlinkは維持する。mainの`.gitmodules`・取得手順・gitlinkと、個人originを正しい保存先へ訂正する。
- mainへ個人macro・ログ・ROOTを追加しない。解析本体・本番出力・未追跡workspacesの内容は変更しない。force push・Git履歴書換えは行わない。

## 検証と残件

統合前HEADは双方の`refs/backup/remote-correction-20261116`で保護した。作業証跡は親の`rootfile/git_remote_correction_20261116_001/`に保存する。統合後のファイル保全、公開境界、remote HEAD、親gitlinkを照合して結果を追記する。

17:14 JST：push前の独立blob監査PASS。旧個人23,354ファイルと正しい保存先2,409ファイルに欠損なし。README・gitignore・当日の記録追記を除き、mode/object IDが同一。`Lambda_check/`・`Nuclear_id/`のtree IDも同一。新規保存23,356ファイルの保存方針検証と既存10 testsがPASS。共有解析の手元差分も統合前と同一である。

既存の正しいmainを基点に統合treeをcommitし、そのtreeを維持して旧個人HEADもmerge祖先として記録する。これにより、既存symlinkを変更せず、通常pushで正しい保存先の履歴と個人作業の履歴を両方保持する。過去資料のCRLF・空白は今回の移行で正規化しない。

旧GitHubリポジトリ自体の削除と、mainから旧参照を除く訂正は別操作である。この環境にはSSHのGit push認証はあるが、GitHub管理API用の認証は確認できない。旧リポジトリ自体の削除は未実施であり、移行完了後にユーザーによるGitHub側の対応を確認する。

tree解析の新規実装計画は、この訂正対応の間は中断する。新しい解析コード・設定の実装、job投入、PicoDst再読は行わない。
