# EP区分数とproton–⁴He原点cutの比較

2021 P24iy・較正388runのrich tree20期間を使う。PicoDst再処理は不要。未回収34producer・未較正81runの保留を維持する。

- [x] ユーザー指定はproton–⁴Heのみ。SE/MEのΛ–⁴He pairにstrict `abs(Δφ*)<0.05 rad && abs(Δη)<0.05` を適用する。event/候補/poolは維持する。
- [x] EP pool N=1/8/16/32を物理Ψ1から分類。保存済みcanonical EP8較正とEP-valid母集団は共通。
- [x] r=1.4m、global momentum η、各元eventの磁場とraw rigidity/chargeで幾何を計算。ME両方向で元eventのcacheを使う。到達不可はkeep+計数、cache欠落/取り違えはエラー。πはQAのみ、³Heは変更しない。
- [x] batch一致build、EP117951 checks、幾何33 checks、実Legacy164339 checksを通過。SE/MEの非空reject、境界、wrap、cache所有、³He不変を確認。
- [x] 同じ40万保存eventの8条件pilot、旧EP8 offとの物理ヒスト一致、cut会計、EP投影和を確認。
- [x] 7新条件×20期間=140jobを投入。EP8 off全統計は既存本番を再利用。depth100・期間内継続・bufferAll・双方向を共通とする。
- [x] 全条件完了後にraw SE/ME scaled・CFを比較。規格化0.4–0.5 GeV/c、背景差引き前。centrality/EPの合成はraw count加算後に規格化。

[第二段階の日本語説明PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/share/figure/20261120_lambda_second_stage/lambda_second_stage_guide_ja.pdf)

[新旧Δφ*–Δη比較PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/share/figure/20261120_lambda_ep_merging_systematics/track_merging_current_vs_old_20261120.pdf)

新tree再読でproton側dip・π側ほぼflatを再現。中心±0.05/周囲の面積補正密度比はproton新0.816±0.021・旧0.821±0.032。bin面積からの診断で、exact cut除外件数/merging効率ではない。

既存EP8 offの0–60%は、規格化後のSE/ME積分比が0.6≤k*<1 GeV/cで約0.9302。現ヒスト上限1 GeV/c。それより上や新条件の改善は未評価。

原本：`rootfile/femto_lambda_ep_merging_systematics_20261120_001/`。source_before_manifest/cases/各test・比較作図helperを保持。今回git commit/pushはしない。

本番7条件140jobを投入・独立検証済み。既存EP8 offを再利用し、全8条件を比較する。完了監視helperの保存障害はユーザーの手動承認への切替回答後に解消し、監視は正常終了（2026-11-21確認）。全件正常終了・ROOT回収・provenance・規格化を確認してCF/plane/raw PDFを作成し、最小メールを1通送る。全統計の比較図は完成。11/21の完了記録を参照。

## 本番投入・検証記録

2026-11-21T00:47:51.825372+09:00 JST。全8pilotで400000event完走。旧EP8 offの1311旧ヒスト・9ROOTの全event台帳が完全一致。許可された7source/doc以外の既存31ファイルもSHA不変。

- ep1_mergingOff: `99BF51C534E3EDB9BF7D4E00D0719DC2`（20job）
- ep1_mergingOn: `A41640F4766ED2F6E154BD81C5343138`（20job）
- ep8_mergingOn: `A730E40D32DC4EDD890D4FFB7F46F217`（20job）
- ep16_mergingOff: `6A21EE258585F2B438551377515209BF`（20job）
- ep16_mergingOn: `EA9E4CBE1B57C76D5337C29F5EA78FF5`（20job）
- ep32_mergingOff: `C44F820BF5267ECB143D2529198B3E3C`（20job）
- ep32_mergingOn: `9F620F840BE41704A1DE1DB1D5082DC4`（20job）

両depth100・期間内mix継続・copyInputLocally・1期間/worker・4096MB RAM/約64GB scratchを全140workerで確認。独立監査：`independent_submission_audit.json` PASS。最長CPU見積31.1h、scheduler見積48h。投入後quota余裕1.209 TB（in-doubt控除後）、追加出力100GBを見込んでも500GB reserve以上。

tmux: `femto_epcut_20261120_001`、aggregate watcher PID 1082664。監視/将来比較出力：`rootfile/femto_lambda_ep_merging_systematics_20261120_001/aggregate_watch_001/`。個別caseの再投入/hadd/自動retryは行わない。完了後の比較PDFは同dirの`cf_comparison/`。

## 2026-11-21完了確認

全140新job正常終了、既存20を含む160ROOT回収、比較PDF3種完成。02:55 EST成功メールと一致。高k*低下はEP32でも残り、cutの高k*窓効果は約−0.07%。[結果と検証](../21/ep_cut_matrix_completion.md)。
