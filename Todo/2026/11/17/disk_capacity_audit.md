# 空き容量確保のための使用量・削除候補調査

**最新状態：ユーザー指示で旧2AF候補ROOT 13,027本を削除済み。約2.073 TB分を削除し、soft余裕は2,249.96 GBへ増えた。比較元PDF/QA ROOTは保持。** [削除結果と再確認Todo](old2AF_retirement_and_merging_recheck.md)。以下は削除前の調査記録を履歴として残す。

削除前の2026-11-17読取り調査。user quota使用量は**2,363.71 GB**、inDoubt差引後soft余裕は**178.21 GB（165.97 GiB）**、hard余裕は677.97 GiB。最大の旧2AF候補保存ROOTだけで**2,072.76 GB**、使用量の87.7%を占める。

この調査では削除・移動・圧縮を行っていない。下記は測定済みの占有量と、実行前に判断/検証が必要な候補。GBは10⁹ bytes、GiBは2³⁰ bytesで区別する。quotaに効く実割当を主に表示し、symlinkは追わず、ROOT directory内のhardlinkを二重計上しない。

## 大きいROOT出力

| 対象 | 実割当GB | 件数 | 判断 |
|---|---:|---:|---|
| 旧2AF候補保存ROOT | 2,072.76 | 13,027 | 旧cache廃止の判断後の削除検討候補。重複ではなく候補payloadは独立。 |
| 現行EP1出力全体 | 60.62 | 4,068 | 確定merged保持。元subjob/中間の削除量は未分離。 |
| 現行EP8出力全体 | 60.62 | 4,067 | 確定merged保持。元subjob/中間の削除量は未分離。 |
| LambdaNuclearId出力 | 26.64 | 8,693 | 既存merge対応と現在参照を確認してから判断。 |

ROOTのdirectory集計は通常ファイル分。directory metadataは別で、NuclearIdはdirectory込み26.678 GB。EP1/EP8の各約9.282 GBの確定mergedも上の各60.6 GBに含まれるため、別途加算しない。

最大対象の正確な範囲：

```text
rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2AF731F2C75544890B4E858354DDE43D_*.root
```

旧2AFは2021 P24iyの候補保存あり本番がquota事故で未完走となったもの。13,027本中2,492本は0 byte、非空は10,535本。予定26,994 jobに対する完全出力ではなく、全候補treeの健全性検査も完了していない。小さいpilotではない。後続1CEはhist-onlyであり、旧候補payloadの代替にはならない。11/11〜12には、元PicoのEvent/Trackと組み合わせてtrack-merging QAを再計算し、10,087成功jobの全96histを統合してΔφ*–Δηの最終図を作成した。原点付近のdipを議論した図の入力であり、未使用ファイルではない。11/13の清掃では明示的に保持していた。現在のEP1/EP8のCF資料に使う新しい確定merged入力とは別である。

**旧candidateStoreを今後使わないと判断する場合の最優先削除候補**。削除後に同じ候補情報が必要になればPicoDstから再生成が必要。新しいtree方式へ移る方針に照らして廃止可否を決める。2AFを廃止する場合、soft余裕は単純計算で約2,250.97 GBになる。

[2AF状態記録](../04/job_status_2AF_20261104.md) / [hist-only再投入](../04/hist_only_relaunch_20261104.md) / [11/13の保持判断](../13/quota_cleanup_and_local_push.md)

[2AFを使用した全統計track-merging作図の記録](../11/track_merging_audit_20261111.md)。候補cacheを再計算へ使った事実と、現在のCF作図入力とは別であることを区別する。作成済みQA ROOT/PDFは元2AFとは別ファイルで保存している。

EP1/EP8は各4,058 jobが終了しqueueは0。ただし最終解析には共通4,023 groupだけを採用した部分haddを使う。採用元ROOTは`subjobs/<jobid>/`へ移動済みで、これは圧縮保存ではない。確定mergedと入力/除外manifestは保持する。元subjobs/中間を削る場合は4,023採用分と各35除外分を区別する必要があり、今回その削減量は未分離。directory全60.6 GBを丸ごと削除候補にしない。

[hadd完了記録](../14/pair_rapidity_hadd_20261114.md)

単体で大きい過去merged ROOTも存在するが、それぞれ異なる設定・統計なので古さだけで不要とは判定しない。

| パス（star-analyzerから相対） | 実割当GB | hardlink数 |
|---|---:|---:|
| `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root` | 14.06 | 1 |
| `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14_merge.root` | 11.02 | 2 |
| `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root` | 10.52 | 1 |

1CEと現在EP1/EP8の確定mergedは各nlink=2。別名hardlink片方の削除だけでは当該データの容量は解放されない。過去2021 partial ROOTは採用PID較正の出典でもあり、参照関係を確認して保持判断する。

## 完了ログと生成jobファイル

| 対象 | 実割当GB | 整理方法/現在の検証状態 |
|---|---:|---|
| Nov13 EP1/EP8 raw stdout/stderr 16,232本 | 4.34 | 圧縮log archive未確認。圧縮・member内容検証後にraw削除 |
| うちhadd採用4,023件×2のログ16,092本 | 4.28 | 最初の整理候補。圧縮archive分は解放量から差引 |
| 除外35件×2のログ140本 | 0.05 | 障害/除外の記録として保持推奨 |
| `job/run/`全体 | 20.45 | 直下88,981本が16.53 GB。下表は内訳であり二重加算しない |

| job/run直下のcsh/list | 実割当GB | 状態 |
|---|---:|---|
| 旧2AF `2AF731F2C75544890B4E858354DDE43D` | 6.88 | runmeta/runtime_bundle/sums_artifactsなし。唯一の実行資料として圧縮保全が先。 |
| NuclearId `523CF569AAF0FFC0623304893ECA9B00` | 4.55 | SUMS tarとrunmetaあり。csh basename/size一致、list全coverage/個々の内容SHAは未検証。 |
| EP1 `FD44A5EB1B0AA22C546F169E86482587` | 2.14 | SUMS/runtime tarあり。全member coverage/内容一致は未検証。 |
| EP8 `DA0DEB4347BD834F1096F8D8798F2133` | 2.14 | SUMS/runtime tarあり。全member coverage/内容一致は未検証。 |

小さいcsh/listでもGPFS上の割当blockが大きく、ファイル内容の合計より大幅に場所を使う。2AFの実行資料6.881 GBは唯一性があるため、ROOTの廃止判断とは別に圧縮保全を先に行う。NuclearId/EP1/EP8は既存SUMS tar等の完全coverage・個々の内容一致を確認してから生ファイルを削除する。今回は削除直前の全検証は未実施。

空になった旧log/err directoryにも約0.33 GBのdirectory blocksが残る。`job/run/`親directory自体は約0.54 GB。空directory削除やdirectory再作成による回収は別の限定作業とし、ここでは実行していない。

Nov13の2 package＋2 zipは合計約0.57 GB。runmeta/runtime/SUMS tar/configlog/joblistlog/入力・除外manifestは再現に必要なので保持。`job/run/runmeta/`全体は約0.56 GB。

11/11・11/13に整理済みの旧ログや84本のtest ROOTは、今回の削減見込みへ再計上しない。新treeの10k検証ROOT/closureは現役の根拠なので保持する。

## その他の領域と容量確保の優先順

- mainのROOT/log/err/job以外は約3.84 GB。source/config/Git/個人記録が中心で、容量確保の主対象にしない。
- git_STAR_3外側のworktree/quarantine/backupsは計約0.45 GB。
- 他の本人PWG領域はgit_STAR_3を除き約6.15 GB。そのうち`/gpfs/mnt/gpfs01/star/pwg/yichikawa/tmp_ana/`が約5.55 GB。内容・現行参照は未評価なので名前だけで削除しない。

1. 旧2AF候補cacheの廃止可否を決める。これが圧倒的に大きい。
2. 完了済み採用jobのログを圧縮保全し、検証後にrawを削除する。
3. archive照合後、生成csh/listを整理する。2AF実行資料は新たな圧縮保全が必要。
4. 必要なら、現行確定mergedを残して、検証済み元subjob/中間だけをmanifest単位で整理する。

新compact treeの外挿は約465.31 GB（433.35 GiB）で、現soft余裕との差は約287.10 GB。ログ＋生成ファイルだけを全て整理してもこの差は埋まらない。旧2AFの廃止、別保存先、または他の大きいROOTの整理が必要になる。外挿は単一run/10k sampleの値で、初段QA・metadata・後段/一時コピー等は別途見込む。

[機械可読の測定値・範囲と制約](../../../../verification/disk_capacity_20261117_001/measurements.json)。この記録は使用量調査の結果であり、削除済みとの報告ではない。追加の網羅的全filesystem scanやGit commit/pushは行っていない。
