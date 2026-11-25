# EP分割・原点cut比較：本番終了確認

2026-11-21（America/New_York）。新規7条件140jobすべて正常終了。既存EP8・cutなし20出力を再利用し、全8条件160ROOTの回収・比較図生成を確認。各条件1,819,447,887eventsで、条件間で共有する母集団。

最終native終了：2026-11-21T02:44:35-05:00。成功通知：2026-11-21T02:55:39.707797-05:00。全160ROOTは91,990,858,528bytes（新140出力80,499,325,024bytes）。独立監査4762checksで160ROOT current stat/header・320log SHA・実完了marker・native入力順・fresh queue0を確認。requested/completed/available一致、Maker completed、chain Warn/EOF/Err/Fatal=0。メール再送なし。

## 結果

0–60%、raw signal・背景差引き前。H=[0.4,0.5)GeV/cで規格化。表はraw窓積分比で、CF点の単純平均ではない。

| EP分割 | C[0.6,1) cutなし | cutあり | C[0.8,1) cutなし | cutあり |
|---:|---:|---:|---:|---:|
| 1 | 0.890103 | 0.889507 | 0.864213 | 0.863588 |
| 8 | 0.930228 | 0.929576 | 0.912158 | 0.911461 |
| 16 | 0.931937 | 0.931301 | 0.914165 | 0.913483 |
| 32 | 0.932518 | 0.931882 | 0.914813 | 0.914130 |

EP1→8で1へ近づくが、16/32細分の追加変化は小さく、全条件で高k*低下が残る。cutによるC[0.6,1)の相対変化は−0.0669〜−0.0701%。この小ささは高k*窓に限定し、最初の0–0.02binでは低統計の約3%変化がある。独立raw再集計3720checks差ゼロ。

SEはcut固定で全EPのnative binが一致。0–60%の除外SE1856/2221644=0.08354%、EP8ME457166/459646635=0.09946%（flow込み）。元軸上限1GeV/c、規格化/共有event/MEの共分散・条件間有意性は未評価。

## 比較PDF

- [CF条件比較](cf_comparison.pdf)
- [各EP binのCF](cf_plane_diagnostic.pdf)
- [scaled SE/ME](raw_se_me_scaled.pdf)

保存先：`share/figure/20261121_lambda_ep_cut_comparison/`。元との同一SHA・全8postprocess artifact SHAを確認し、0–60%のCF/raw代表pageを目視。postprocessは160入力・226曲線、input/pins unchanged、treeEntriesRead0。

出力回収後のsoft quota余裕は1.135TB（in-doubt控除後）、graceなし。削除/hadd/再投入/メール再送なし。初段34未回収/未較正81runの保留は維持。

原本：`rootfile/femto_lambda_ep_merging_systematics_20261120_001/aggregate_watch_001/`、今回監査：同studyの`completion_audit_20261121/`。
