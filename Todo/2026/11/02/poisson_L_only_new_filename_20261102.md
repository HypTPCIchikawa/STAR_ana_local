# L結果の総合レポートを新しい名前で出力

状態：完了。ユーザーの閲覧では元リンクが34頁と表示されているため、元パスには今回一切書き込まず、L結果だけの20頁レポートを別名で出力した。ローカルでは元パスも20頁であり、表示の違いの原因は未特定。

- [新しい総合PDF（L only、20頁）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/summary_4He_KFParticle_3p85_PoissonL_only_20261102_v2.pdf)
- [新規出力・元パス保持の検証記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/summary_4He_KFParticle_3p85_PoissonL_only_20261102_v2.json)

通常fitと非負定数BG再fitはともに`SLQNR0`。非L結果の図は含めない。全20頁・15図頁のL表記と元図のhashは先行検証で確認済み。独立したページ数確認も20頁。

Case2は新レポート13頁。SE method3、50–60%、k*=0.26–0.28 GeV/c。BGは定数`0.14694668995939236 counts/bin`（傾き・二次係数0）。fit範囲の70bin中57binが0、合計27countsである。

元のPDFパスについては作業前後のSHA256・byte数・mtimeNsがすべて一致することを確認した。数値fitの変更、production/hadd、Git commit/pushは行っていない。

独立したPDF文字抽出でも、保存された旧34頁版の27頁が旧χ²のCase2に対応し、新20頁版の13頁がLのCase2に対応することを確認した。ユーザーが見た27頁は旧版の特徴と一致する。表示経路やキャッシュが原因かは断定していない。

追記：高k*の低下と規格化の解釈は[ME/k*監査](cf_mixing_kstar_audit_20261102.md)を参照。0.4–0.5は比較用スケールであり、物理的plateauとして認定していない。保存された数値・PDFは変更しない。
