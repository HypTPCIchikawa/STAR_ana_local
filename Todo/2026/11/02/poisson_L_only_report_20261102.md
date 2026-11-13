# L付きfit結果だけの総合レポート

状態：完了。直前24頁PDFにχ²/L比較図3件（4頁）と旧値の列が残っていたため、L結果だけの構成に改める。

- [x] group/fineのfitOptionsがともにSLQNR0であることを確認。
- [x] 非L結果を含む3図を特定。
- [x] 旧比較図・旧純度列・toyの旧χ²値を除外。
- [x] 全図にPoisson L (SLQNR0)を明記。
- [x] 15図頁を含む総合PDFを再生成し全頁を確認。
- [x] 指定PDFを差し替え、直前版を保存。
- [x] README・日報・個人snapshotを更新。

fitの再実行や数値変更は行わず、採用済み004のL結果を使う。Method3/5の9中央度純度、Method3/5/6の2・3分割CFとSE/ME純度、代表9例を残す。

## 現在の総合レポート

- [L結果だけの総合PDF（20頁、新しい名前）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/summary_4He_KFParticle_3p85_PoissonL_only_20261102_v2.pdf)
- [直前の比較を含む24頁版（保存）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_pid_tof_updated_20261030_001/pdf/revisions/before_L_only_20261102_001/summary_4He_KFParticle_3p85.pdf)
- [設定・全図由来・検証](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_poisson_likelihood_20261102_004/pdf/l_only_report_20261102_001/README.md)

全20頁のヘッダーと全15図頁のキャプションにPoisson L (SLQNR0)を明記した。本文の旧P(SE)列・toyのχ²推定値も除外した。Method3/5の9中央度純度、Method3/5/6の2分割・3分割CFとSE/ME純度、代表9例のmass fitを掲載する。全20頁を画像化して目視確認、非L比較図ゼロ、LaTeXはみ出し警告ゼロ。独立した図・内容監査もPASS。

元のL数値、質量分布、fit結果、純度/CF TSV、規格化表は変更していない。004の1960点・632定数refit・747有効CF点の検証を引き継ぎ、レポートへの採用図と表記だけを改めた。現在のPDFのSHA256は`a1e1e2410ecef0e4fc7d4057d3c378c88b7ac2fef71a863fde7f2e05956f2bde`。

[最終照合・個人snapshot保存記録](L_only_final_checks_20261102.json)。小さな成果物と生成コードを個人管理へ保存し、全コピーのhashを照合した。大きなROOT・数値配列・ビルド生成物は既存policyで除外し、元ファイルは保持。commit/pushは行っていない。

ユーザーの表示相違への対応として、[新しい名前へ出力](poisson_L_only_new_filename_20261102.md)。元のPDFパスは今回変更していない。
