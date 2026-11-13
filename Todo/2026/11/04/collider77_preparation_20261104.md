# 7.7 GeV collider解析の非干渉準備

最新状態：2026-11-06のユーザー指示で準備を再開。[本日のTodo](../06/collider77_preparation_20261106.md)で進捗と検証結果を管理する。

以下の保留は11月4日時点の履歴：現行3.85 GeV jobの終了を確認してから再開する方針とし、workspace作成・変更・build・投入は行わなかった。11月6日はjob終了済みだがhaddが稼働しているため、独立workspaceで準備を進める。

非干渉で準備する方法は確認したが、実行用workspaceの作成・7.7 GeV用設定の検証・job投入は未実施。対象はユーザー指定の2021年P22ibだが、catalogでの実在・正式なtrgsetupnameは未確認。

## 現行jobとの分離

実行中の3.85 GeV解析はSUMS ID `1CE33994716F7896423A80AB994BAB14`。workerが起動するたびに、元の `analysis/config/lib/include/StMaker/StRoot` をscratchへコピーする。調査時点でも全件生成・投入を継続中であり、元のruntimeを変更してはならない。

別の物理workspaceなら準備可能。例として、現在の `star-analyzer` と同階層の `star-analyzer_auau7p7_2021_P22ib` を使用する。これは予定名であり、まだ作成していない。

- 通常コピーで必要なソース・設定を分離し、hardlinkや元runtimeへ向くsymlinkは共有しない。Git未追跡の設定・個人ソースも必要なので、Git worktreeだけで完結したと考えない。
- `lib`、yaml-cpp build、`.current_mainconf`、ACLiC用tmp、job/run、log、err、rootfileを全て独立させる。
- コピー済みlibraryに元の `lib` への絶対RPATHが残るため、別workspaceでbatch-matched toolchainにより再buildする。継承したlibrary検索pathも確認する。
- `analysis_info` の `anaName/workDir/logDir/errDir` を別解析用に揃える。joblist生成・submitは別workspaceだけで行う。
- rootfile・成功ログ・既存job成果物をworkspaceへ大量コピーしない。少統計検証前に容量を再評価する。

読み取り監査で確認した共有更新箇所：`script/setup.sh:68`（current mainconf）、`Makefile:105–106,241–249,344–347`（library/yaml-cpp build）、`analysis/run_anaFemtoLambda_4He.C:55–60`（ACLiC tmp）、`script/analysis_info_helper.py:623–631`（joblist）、`job/run/setup_symlinks.sh:9–12`（runtime link）。絶対RPATHの設定は `Makefile:308,327,335`。

runtime・script・template・testsの候補ファイルは約23.6 MB、GPFS実割当約245 MB。yaml-cpp source/buildは未計測であり、準備全体や解析出力の容量見積もりではない。

## 準備Todo

- [x] 現行jobのlive runtime依存と共有書き込み先を調べる。
- [x] 元workspaceでbuild/setup/joblist生成を行わず、独立workspaceを使う方針を整理する。
- [ ] 2021 P22ibの正式なdataset・library・PicoDst schemaをcatalogで確認する。
- [ ] 必要な未追跡ファイルも含め、独立workspaceを作成してpath/linkを監査する。
- [ ] 7.7 GeV用mainconf・analysis_infoを用意する。8次元 `hNuclearPidPairJoint` は無効化を維持する。
- [ ] collider用の頂点Vz/Vr・trigger・run/year別centralityを確認する。固定標的用の頂点付近やMEのVz clampを流用しない。
- [ ] KFのPV幾何・固定標的の仮定、lab/CMS rapidityのshift、ヒストのVz/eta/Y範囲を確認する。
- [ ] 4HeのdE/dx/TOF PIDをP22ibの実データで点検する。3.85 GeVで調整した係数の妥当性を前提にしない。
- [ ] 独立領域でbatch-matched buildと少統計実行を行い、event/pair/histのfillと8D不在を確認する。
- [ ] 少統計の出力量とcatalog統計から容量を見積もり、本番投入条件をまとめる。

既存の13.5 GeV比較候補は `config/mainconf/main_auau13p5_anaFemtoLambda_4He_KFParticle.yaml` および対応するevent/centrality/KF設定。7.7 GeV専用設定は今回のファイル一覧では見つからなかった。13.5 GeV設定の中身の詳細確認は未完了であり、そのまま使えるとは判定していない。

## 未実行のcatalog確認

```bash
get_file_list.pl -nocache \
  -keys trgsetupname,production,year,library,filetype,path,filename,events,size \
  -cond production=P22ib,trgsetupname~7p7,year=2021,filetype=daq_reco_picoDst,storage=local \
  -limit 10
```

helpで上記option/keyは確認したが、catalog照会本体は自動承認レビューのタイムアウトでdispatchされなかった。設定読取も同じ理由で一部未実行、通常sandboxの読取はnamespace作成エラーで開始できなかった。これらをデータなし・設定不備とは扱わない。未実行の照会結果や妥当な頂点cut値は記録していない。

今回の操作は読み取り調査とこの個人Todoの作成のみ。現行3.85 GeVのruntime/config/libraryとjob運用は変更せず、7.7 GeVの本番投入も行っていない。
