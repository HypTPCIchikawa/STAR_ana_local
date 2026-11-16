# Pair rapidity 10区分・EP対象run・quota確認

**最新状態：22:31 JST、EP1/EP8各4,058jobs投入完了、完了メール監視起動済み。81runは保留。過去の未投入記録は作業履歴として保持。**

- [x] ユーザー方針：未較正81runの解析と較正・共有補正試験は保留。較正済388runでEP分割の有用性を確認できたら81runの較正を検討する。現時点で81runの処理を再開しない。
- [x] 物理Lambda＋He4のlab pair rapidityを−1.5〜0.5の10区分でSE/ME k*へ保存。有限範囲外はflow、既存He4単体YLab分割は維持。
- [x] 8本TH3Dでfull/signal/leftSB/rightSB×SE/MEを保持し、隣接Y区分とcent9を後処理でmerge可能にする。
- [x] 生成YAML・compiled Maker・QA・意味ある回帰・実Picoで既存SE/ME不変と新区分閉包を検証。
- [x] 最新user/groupのblock/file quotaと本番出力/ログ/一時hadd増分を見積もる。1 EP8campaignとEP1/EP8比較2campaignを区別。

本番jobの投入はまだ実施していない。設定・quotaの確認は投入後の容量を見越した事前監査として行う。


## 実装・検証完了

[詳細と再現記録](../../../../analysisnote/20261113/pair_rapidity_implementation_and_quota.md)。約160万件の回帰(1603559)、projection130、実Pico10000、旧1707physics hist完全一致。25files/jobなら2系列+haddの予測166.9GB/stress219.2GBが現soft余裕286.1GB以内。388run73054files・約2923jobs/系列の専用XML準備、未投入。

- [x] 実worker完了後のruntime/staging・実ROOT容量・quota・共通group確認と本番比較を完了。共通4023groupの確定haddと、[2026-11-15のEP/rapidity比較PDF](../../../../analysisnote/20261115/ep_rapidity_highk_cf_report.md)を参照。初期欠損groupの問題は下記の履歴として保持。
- [ ] EP分割の有用性を確認できたら、保留中81runの較正を検討する。現時点は保留維持。

- [x] 21:40 JST 再開後の最終QA目視・独立review・最終差分/SHA/manifestを完了。quota再確認でも2系列+stressがsoft内。本番job未投入、81run保留。


## 22:22 JST：本番投入を開始（過去の「未投入」記録を更新）

- [x] ユーザー依頼により、較正済388run・73,054filesのEP1／EP8本番を22:05 JSTに開始。未較正81runは保留維持。
- [x] SUMSの実分割4,058jobs/系列見込みでquotaを再計算。2系列＋hadd込みstress256.35GB、開始時soft余裕286.00GB内。
- [x] 切断後も投入を続け、runmeta保存後に最小内容の完了メールを各系列1通送る永続ランチャーを起動。
- [x] 全jobの投入完了、EP1/8実入力group一致、メールwatcherの実起動を確認する。
- [x] 回収後の実ROOTサイズ・live quotaをhadd復旧時に確認。EP1/8 CF・rapidity比較と日本語PDFを2026-11-15に完了（上記リンク）。

[投入の詳細・SUMS ID・状態](../../../../analysisnote/20261113/pair_rapidity_production.md)。


## 22:33 JST：投入・監視確認完了、解析完了待ち

- [x] 各系列4,058jobs、73,054filesのcoverageと全group順序一致を確認。各系列のrunmeta/configlogを保存。
- [x] 完了メールwatcher2本を実起動し、各4058jobs・宛先・minimal本文・process生存を確認。
- [ ] 初期EP1の5workersに見つかったstaged PicoのZombieエラーを回収後に調べ、同じ入力groupで回復する。外側Condor exit0だけでは成功判定しない。
- [x] 完了ROOT/ログを確認し、両系列の共通成功4023groupに揃えてhadd・CF比較を完了。未出力35groupを含む全投入統計とはしない。

[最終状態・初期エラー詳細](../../../../analysisnote/20261113/pair_rapidity_production.md)。現在は解析実行中であり、全統計取得済みとはしない。


## 2026-11-14 10:59 JST：job終端・共通haddへ

両4058jobs終了。ROOT未出力はEP8 35/EP1 34で、共通成功4023groupの検証付きhaddをtmuxで開始。[最新Todo](../14/pair_rapidity_hadd_20261114.md)。EP1のみ有効3496は元ROOTを保持。全投入統計とはしない。81run保留を継続。
