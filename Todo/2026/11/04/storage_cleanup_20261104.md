# ROOT・成功ログの容量整理（2026-11-04）

最新：13.5 GeV分も明示承認後に整理完了。累計 **419.177 GB**を解放、soft quotaまでの空きは約 **373.356 GB**。下の初回3.85 GeV表は当時の値、追加分は末尾参照。

初回3.85 GeV分の完了時刻：2026-11-04T13:33:11.760707+09:00。対象は3.85 GeVの検証済み旧解析04BF/0BEF。ユーザーの指示により、統合済み入力ROOT・hadd中間ROOT・成功済みstdoutとSUMS生成物を削除した。

## 実測結果

| 対象 | 削除数 | 論理容量GB | 解放した割当容量GB |
|---|---:|---:|---:|
| 統合済みsubjob ROOT | 53,911 | 289.856 | 297.212 |
| hadd中間ROOT | 552 | 103.050 | 77.681 |
| 成功stdout・job生成物 | 162,216 | 3.781 | 42.532 |
| 合計 | 216,679 | 396.688 | **417.425** |

GB=10^9 bytes。`mmlsquota -Y`の使用量差も**417.425GB**で、削除時に記録した最後のlinkのst_blocks合計と完全一致。GPFSの小さいファイルの最小割当が262,144bytesのため、ログ/jobは論理容量より大きく解放された。中間ROOTの5つは保持した最終ROOTへのhardlinkであり、その別名削除を容量解放には数えていない。

- quota使用量：2.643TB → **2.226TB**。
- soft quota：2.597TB、hard limit：3.147TB。
- soft quotaまでの空き：**371.622GB**、hard limitまで：921.378GB。猶予期限切れの超過状態は解消した。

## 保持した成果物と削除の根拠

04BFは26,955入力を統合した通常1,786ヒストの最終ROOTと、3つの8D Sparse sidecarを保持した。0BEFは26,956入力の採用済み部分統合ROOTを保持した。全5ROOTはaccepted fingerprintと一致し、ROOT5で非Zombie・非Recovered、true/mix内部を含む再帰ヒスト数とSparse/split manifestを確認。0BEFは全ファイルSHA256も採用済み検証値と一致した。

削除入力はarchive_manifestのsize/inode/device/mtimeと実行直前に照合した。成功log/jobは最終mergeのselected indexに限定し、csh/listを保存済みSUMS圧縮archiveへ全件収録済みと確認した。失敗・除外job、設定・runmeta・session/report・検証記録・小さな圧縮provenanceは保持した。errディレクトリは今回の削除対象外。

今回2AFの未統合ROOT（13,027本、約2.071TB）と失敗調査・再投入に必要な記録は削除していない。旧auau3p85_anaLambdaNuclearIdの8,624本は、既存mergeとの入力対応を確認できず保持。13.5 GeV分は当初の承認待ちを経て、下記の明示承認後の追加整理を完了した。

## 今回jobの全件保存見込み

正常完了ログとROOTヘッダーfEND一致を満たす21本を使用。入力event数とROOTサイズの関係は次のとおり。

`ROOT bytes/job = 4,562,918.83 + 2,449.780609 × 実入力event数`

R²=0.999538、残差RMSE=0.997MB。catalog記載2,184,032,750 eventsに、標本の実入力/catalog比0.96404959を適用すると、26,994jobsのROOTは**5.281TB**。同じ2021 P24iyの旧解析の実測入力2,100,726,543 eventsで計算しても5.269TB。

| 保存内容・条件 | 見込み |
|---|---:|
| 今回のROOT全件 | **約5.3TB（目安5.2–5.5TB）** |
| 現存約2.071TBからの単純差 | 約3.210TB |
| 全件を新jobidで再実行し、現存分も保持する場合 | 約7.353TB |
| log＋job生成物の割当量 | 約21.300GB＋小さな再現用archive等 |
| stderrも含めた割当量 | 約28.376GB＋再現用archive等 |

目安幅は信頼区間ではなく、正常標本とevent数補正に基づく容量計画用の幅。全ROOTの健全性は未検査であり、単純差は既存出力をそのまま有効と認定した値ではない。未処理runでの候補数分布の変化、再投入で生じる重複保存やhaddの一時領域は別途考える。

実ROOT3本では候補保存TTreeの圧縮bytesが全体の92.3–97.1%。今回の大容量の主因は8Dヒストではなく、物理候補とpre-PID QA候補を保存するTTreeである。今の保存設定で全件を置くには、整理後のsoft quota空き371.622GBに対してなお約2.838TB不足する。解析コード・保存設定・PIDをこの整理で変更しておらず、再投入も行っていない。

## 記録

- [削除・容量のsummary](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/cleanup_summary.json)
- [容量モデル・quota実測](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/storage_estimate.json)
- [保持ROOT検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/retained_root_check.json)
- [候補TTree容量の実測](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/candidate_storage_components.json)
- [全削除manifest・実行journal・検証資料の圧縮保存](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/cleanup_audit_records.tar.gz)

最初の補助検査は最上位のみを数えて145ヒストとしていたため不一致となったが、true/mixディレクトリを再帰的に数えると1,786で一致した。データ異常ではなく検査範囲の修正。容量component JSONのdenseHistogramsは最上位のみの補助値で、TTree比率と全ヒスト数の検証は区別する。

## 13.5 GeV分の追加整理（明示承認後）

承認：ユーザー回答「13.5 GeV分も追加で削除」。完了：2026-11-04T13:39:59.084342+09:00。対象SUMS ID `1213058CE62BCF113B2E3C75ACBD7053`。

| 対象 | 削除数 | 論理GB | 解放割当GB |
|---|---:|---:|---:|
| 統合済み元ROOT | 542 | 0.982 | 1.068 |
| hadd中間ROOT | 7 | 0.498 | 0.257 |
| 成功ログ・生成物 | 1,630 | 0.036 | 0.427 |
| 追加合計 | 2,179 | 1.516 | **1.752** |

13.5 GeVの統合ROOTは削除前後ともSHA256 `50c9b4b3d7742b0d654c663eab81901f10119c340daa17993d7f740bb99da03b` で、採用済み検証と一致。元ROOT542本と中間7本を削除し、最終ROOTへのhardlink別名は解放量に重複計上していない。除外149・250・251のcsh/list全6本、失敗時作業記録、SUMS圧縮archiveを保持した。成功csh/listの全収録とgzip CRCを確認してから削除。独立検査もPASS。

追加のquota実測差は1.746141GB。累計は削除時の割当記録で**419.177GB**（218,858ファイル）。小さな監査資料・レポートの追加保存があるため、最初から現在までのquota純減とはその分だけ異なる。今回jobのROOT全件約5.3TBという見積もりは変わらない。

[追加削除・検証記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/completed.json)／[独立検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/independent_postcheck.json)／[累計summary](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/storage_cleanup_20261104_001/cumulative_summary.json)。
