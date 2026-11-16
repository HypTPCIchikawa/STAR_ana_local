# 新1CEによる Λ–⁴He 高k*系統性（2026-11-11 JST）

新試料でも高k*低下は再現した。nσ窓を狭める効果は中央試料で小さく、rigidity・mT依存が大きい。原因の特定、最終CF規格化、PIDの新規採用は行っていない。

[結果と図付き報告](lambda_he4_highk_systematics_20261111.pdf)／[18頁atlas](cf_highk_systematics_20261111.pdf)／[全条件の数値表](metrics.tsv)／[計算・設定](results.json)／[独立検証](independent_validation.json)。

値はk*=0.7–0.8 GeV/cの積分raw CF、0.4–0.5 GeV/cで比較用規格化。純度補正なし。誤差は保存countに基づく統計誤差。

|条件|0–20%|20–80%|
|---|---:|---:|
|全選択⁴He|0.9079 ± 0.0033|0.8942 ± 0.0083|
|nσ ±0.6|0.9093 ± 0.0049|0.9161 ± 0.0128|
|nσ ±1.0|0.9065 ± 0.0040|0.9011 ± 0.0102|
|R<1.0|0.8576 ± 0.0051|0.8816 ± 0.0101|
|1.0≤R<1.5|0.9152 ± 0.0057|0.9128 ± 0.0170|
|1.5≤R<2.0|0.9699 ± 0.0071|0.9822 ± 0.0333|
|2.0≤R<5.0|1.0002 ± 0.0165|1.0545 ± 0.1320|
|2.4≤mT<2.6|0.8899 ± 0.0052|0.8929 ± 0.0107|
|2.8≤mT<3.0|0.9667 ± 0.0106|1.0277 ± 0.0559|

0–20%のnσ±1とinclusiveの差は−0.00139±0.00222（nested overlapを考慮）。R≥2はTOF+m²窓を要求するため、R依存をTOF単独の原因としない。Rに相関する混入も、この独立1D投影だけでは除外しない。

³He物理質量／⁴He質量仮説の0–20%値は0.95828／0.95777で、⁴Heの全低下を質量変更だけでは再現しない。上下Λ信号窓は同様に低下し、左右sidebandはより平坦。旧04BFと新1CEの0–20%値は0.90781／0.90788。新旧はイベントが強く重なるため独立統計の差として有意性を付けない。

## 入力・検証

2021 P24iy、SUMS `1CE33994716F7896423A80AB994BAB14`。mainconf `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`。runmeta/configlog/merge resultは同IDで結合する。26,945/26,994、49件除外の部分統合。入力2,099,784,951、選択2,066,742,256 events。

新旧1,656histをROOT5で抽出し入力不変を検証。1,680系列、84,000表示bin、1,176投影・mass半窓検査。独立に318数値と24導関数を照合。³He ME overflowのみTH2F丸め差最大31count（2.52e−7相対）、0≤k*<1の結果へ影響なし。通常binの検査許容差は維持し、flowだけFloat ULPで別検査した。

[主検証](validation.json)／[新抽出検証](extraction_validation.json)／[旧抽出検証](extraction_validation_04BF.json)／[Float独立監査](legacy_float_independent_audit.json)／[描画記録](plot.log)。

## 再現

数値計算はホストのシステムPython/NumPyを使用する。CVMFS版NumPyはホストでlibgfortran不足なので、PYTHONPATHを除いて既存システム環境を使う。

```bash
env -u PYTHONPATH OPENBLAS_NUM_THREADS=1 /usr/bin/python3 STUDY/source/analyze_highk.py --settings STUDY/settings.json
```

抽出コマンドは`extraction_command.json`と`extraction_command_04BF.json`、入力指紋は`extraction_provenance*.json`に保存。描画は同mainconfを指定した`script/singularity_root_batch.sh ... root4star -l -b -q 'STUDY/source/DrawHighK.C+("STUDY/points.tsv","NEW_ATLAS.pdf")'`。新しい出力名を指定する。設定値は`settings.json`にあり、production設定を変更しない。

## 解釈と残作業

規格化0.4–0.5は物理plateauとして認定しない。count由来の規格化共分散とnested cut相関を含め、event/shared ME pair相関・純度やmassモデル系統は含めない。傾きは0.5–0.8をfull covarianceで記述し、低統計の欠測binは採用せず実際のfit binをmetricsに記録する。

candidate cacheとjoint8D Sparseがないため、R×nσの同時条件、低R TOF別CF、multiplicity/Vz/run/磁場/EPの再混合は未実施。prePIDは元physical poolのQAであり、quality/R/TOFを残す。³He core→⁴He混入率の定量や最終純度補正、controlled remixは未完了。

次の切り分けは、R/mT分布を合わせた比較と、候補/PIDを固定したcontrolled remix。後者には保存候補または別途PicoDst再処理が必要。今回、追加production・メール・git commit/pushは行っていない。

atlas頁: 1–2 species、3–4 mass窓、5–7 nσ、8–10 R、11–12 mT、13–14 prePID、15–16 normalization、17–18 production。先頭2頁を加えた報告PDFではatlas頁に2を足す。
