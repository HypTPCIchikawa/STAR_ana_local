# Tree本番のrun時期別haddと未回収34件

**最新：2026-11-19 00:58:44 JSTに全20群のhadd・出力照合・全体coverage検証を完了。** 正式ROOT20本、合計465,104,188,956 bytes、`.partial.root`残存0。4,024 producer・保存1,819,447,887 eventsを保持。原本は保持し、第二段階SE/ME解析は未実行。下の起動時の記録は履歴として保持する。

ユーザー依頼：取得時期の近いrunを10〜30個のROOTにまとめ、後段SE/ME解析の入力を準備する。元producerを任意のevent位置で分断せず、run番号順の近接区間を20群にする案を検証する。

## 初段jobの最終状態

2026-11-18 17:32 JSTの通知記録で全4,058 jobsがterminal、queue0。最後のsubjob2869もrequested=completed550,829、reconstructed543,693、status0。回収は**4,024 ROOT**、未回収は従来と同じ34件。メールは17:39 JSTにlocal sendmail受付済み。SUMSの外側exit0だけでは解析正常終了の証拠にしない。

- mainconf：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml`
- SUMS：`9037A6F05FAD08FC8E9BD9286C8A38A5`
- runmeta：`job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_tree_treeWrite20261117_9037A6F05FAD08FC8E9BD9286C8A38A5.json`
- configlog・入力・runtime保存：`rootfile/femto_lambda_tree_production_20261117_001/`と同runmeta参照先。
- 統合作業証跡：`rootfile/femto_lambda_tree_grouped_merge_20261118_001/`。
- final snapshot：正常完了ログからinput1,842,806,823、保存event1,819,447,887。ROOT論理471,489,115,031 bytes、割当472,021,729,280 bytes。これは未回収34件を除いた統計。

## 検証・残作業

- [x] 最後のjobを含む4,024原本を固定し、全成功footerとfile statを照合。最初の4,023本のsnapshotも履歴として保持。
- [x] ROOT5/SL24yで辞書、fast merge、圧縮6保持、同process内MaxTreeSize1TBを確認。
- [x] 2source実merge→check：全tree件数/schema・metadata値/順、108QA全bin/content/error/entries一致。原本不変。
- [x] 直前quota空き1,773,689,372,672 bytes（17:59 JST、inDoubt控除後）。原本保持の追加統合を予算化。
- [x] 10/20/30案を比較し20群を採用。18:09 JSTにtmux `he4_tree_period_hadd_20261118_001`で逐次統合driverを起動。各群185〜240原本、入力容量22.805〜24.236 GB。開始時は最初の群のsource catalog検証中。
- [x] 各群catalog・統合・独立照合・原本statの確認後に正式名へ確定。2026-11-19 00:58:44 JSTに全20群完了、全体coverage・4,024 unique producer segment・原本stat不変を確認。各群108QAの全bin/content/error/Entriesの照合を通過。
- [x] 未回収34件は前回EP1/EP8と入力group・エラー分類一致、33件は同一問題Picoを確認。612Pico中309 positive＋1 zeroで実測合計8,594,396 input entries、302未知。欠損入力割合の下限0.464210%、run別file欠損最大36.69%。accepted/pair欠損は未確定で、無視できるとは断定しない。[監査資料](../../../../analysisnote/20261118/tree_missing_inputs_20261118.md)。
- [ ] **保留（2026-11-19ユーザー指示）**：[未回収34件の原因解明](../19/unrecovered34_cause_investigation.md)。原因の切り分けと再処理要否の判断は将来課題とし、現在は追加調査・再投入を行わない。その他579本は正常入力と証明された意味ではない。
- [x] [期間内連続ME・pool100専用設定](../19/period_mix100_production.md)を用意し、40万eventで検証後、11/19 16:10 JSTに20job投入確認。専用設定はfalse/falseと両depth100。元replay標準true/trueとdepth20は再現用に保持する。
- [ ] [Δφ*–Δηの全統計再作成と旧PDFとの傾向照合](../17/old2AF_retirement_and_merging_recheck.md)。今回のIO試験で代替しない。

統合前の原本・native入力list・runtime/較正記録は保持する。未較正81runは保留のまま。後段jobはこのhadd依頼では投入していない。

## 採用分割・進捗の参照先

[20群のrun・容量表と比較](../../../../analysisnote/20261118/tree_run_period_grouping_20261118.md)。元producerを一切分割しないため18runに隣接群へのbridgeが残る。minority側175Picoは回収72,442Picoの0.2416%。run番号順で取得時期を近づけ、ROOT生成mtimeを取得時刻に使わない。

- 本番driver設定：`rootfile/femto_lambda_tree_grouped_merge_20261118_001/merge_config.json`。1並列、fast/no local copies、圧縮6、同ROOT processでMaxTreeSize1TB、float QA許容0 ULP。
- live状態：同studyの`merge_status.json`、`tmux_merge.log`、`operations/group_XX/{catalog,merge,check,published}.json`。
- 起動時の履歴：最初のgroup00は236sourceのcatalogと独立snapshot照合を通過し、18:12:55 JSTに実merge開始。この時点では未完了だったが、その後上記時刻に全群を完了。
- 正式出力：`rootfile/auau3p85_anaFemtoLambda_4He_tree/run_periods20261118_001/`。検証中は`.partial.root`、検証後だけ正式名へ確定。全群完了後に`grouped_roots.list`と`completed.json`を保存する。
- 群ごとのmetadata/all tree件数/schema・QA全bin/content/error/entries、source statを照合。event payload全件の物理再計算ではなく、各segment端点3点と公式Reader少数sampleによるIO検証である。

[個人submodule内の起動・検証証跡](../../../../verification/tree_run_period_hadd_20261118_001/README.md)。共有mainへの解析マクロ・bulk logのcommit/pushはしていない。

## 完了状態の確認

2026-11-19 14:42 JST照会で`merge_status.json`/`completed.json`がともに`phase:complete`、正式ROOT20本、partial0、`grouped_roots.list`の存在を確認。tmuxは処理終了により消滅している。`all_sources_validation.json`はPASS、input1,842,806,823・accepted1,819,447,887・source/segment4,024、原本fingerprint不変。検証scopeは上記のbounded IO確認で、全eventの物理再計算ではない。完了通知メールは今回のgrouped haddには設定していない。
