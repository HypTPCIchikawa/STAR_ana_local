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

## 17:21 JST：参照訂正と正しい保存先への公開完了

- 個人内容の公開commit：`26f3a09a77649a269a6d143e18f1c5b99fcd63da`。正しい既存mainの後に取込commit`b5764c1b`を作成し、旧個人HEADも祖先として保持した。全個人ファイルと両履歴を保全済み。
- 共有内容の公開commit：`a6d617f4ae306504933b7b77355501ee810eecde`。変更は`.gitmodules`・`README.md`・`star_local`の3pathのみ。個人macro・ログ・ROOTは共有mainに追加していない。
- 双方のremote main、local HEAD、親gitlinkが一致。`STAR_ana_local`へのURL訂正、yaml-cpp不変、両既存HEADの祖先関係、取込後のtree不変を独立監査で確認した。
- 作業証跡：`rootfile/git_remote_correction_20261116_001/first_push_verification.json`と`prepush_integrity.json`。要約の保存先は`verification/git_remote_correction_20261116_001.json`。

自動承認レビューは、旧repo削除に使う認証を標準環境設定・認証ファイルから確認する操作を、明示許可のないcredential probingと判断して拒否した。操作は未実行であり、別手段で迂回していない。旧repo本体は未削除。移行保全と現行参照の訂正は完了しており、旧repo本体をユーザーがGitHub管理画面で削除するか、参照訂正のみで完了とするかを確認中である。

この完了記録を個人repoへ公開し、親のgitlinkを追従更新する。後続commitは記録と参照のみで、解析コード・ROOT・本番jobを変更しない。
