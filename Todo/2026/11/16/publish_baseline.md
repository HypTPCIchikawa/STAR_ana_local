# 11/13本番ベースラインの公開

ユーザー方針：EP8各planeのk*保存は実装せず、11/13本番投入時の機能を維持する。未実装draftを退避し、mainと個人サブモジュールをpushする。生成ROOTの後処理macro・個別設定・報告資産・不要ログをmainへcommitしない。

- [x] ユーザーに戻す範囲を確認。C++本体は未変更であることを確認。
- [x] 投入packageとのStMaker/include/analysis 160ファイル、既存config411ファイルのSHA一致を監査。160には公開対象外の生成.so/.dも含む。
- [x] 未実装draft9ファイルを[drafts](../../../../drafts/ep_kstar_20261116_001/README.md)へSHA付きで退避。本番設定・解析本体は維持。
- [x] STAR SL24y/ROOT5 batch整合clean build、共通job管理66 tests、skills同期チェックPASS。
- [x] main公開候補33ファイルを明示し、ROOT後処理macro・tools・個別config・ログ・ROOT/so等を除外。private側605pathを検証してsnapshot更新、535fileを新規/更新コピー。
- [x] サブモジュールとmainをcommit/pushし、remoteとgitlinkの一致を確認。内容公開commitはmain `00bf4bd` / private `37d798a7`（15:59 JST照合PASS）。

[解析記録](../../../../analysisnote/20261116/publish_baseline.md)。新規本番job・Pico再読・メールは行わない。未較正81runは保留。

16:00 JST完了記録を保存。公開境界の独立監査PASS、私用verification10 testsもPASS。検証記録を格納する後続commitと親gitlink更新を加えるが、解析本体の変更は加えない。
