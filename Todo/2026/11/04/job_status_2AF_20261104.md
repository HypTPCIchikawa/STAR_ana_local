# STAR解析jobの状態（2026-11-04、日本時間）

SUMS ID: 2AF731F2C75544890B4E858354DDE43D

確認時刻: 2026-11-04T12:27:38.339367+09:00

キューにjobは残っていないが、解析は正常完了していない。予定26,994件に対し、Condorの終了記録は14,262件（exit0:11,446、exit1:2,816）。ROOTファイルは13,027件、欠落13,967件。全ROOTのヘッダー・内容を検査済みとはしない。

GPFSのuser soft quotaを超えて猶予期限が切れている。後半のworker0.csh/最後のCondor groupは384byteで切れており、reportは1byte、runmetaは未作成。投入後の記録保存と通知worker起動もDisk quota exceededで失敗した。キューが空であることを全件正常完了としない。

容量確保後、正常出力の保持と未完了入力の特定・再投入が必要。再投入・hadd・解析データ削除は行っていない。

状態記録: /tmp/yichikawa/star_job_status_2AF731F2C75544890B4E858354DDE43D_20261104/status.json

最終出力監査：非空・root magic確認10,535件、空または不正ヘッダー2,492件。正常完了ログmarker10,223件。全ROOT内容・候補treeは未検査。詳細を含むメール送信は自動承認で実行前に拒否された。
