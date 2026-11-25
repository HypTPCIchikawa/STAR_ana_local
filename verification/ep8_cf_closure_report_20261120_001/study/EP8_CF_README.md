# EP8・cut Off の CF 図と数値

第一段階treeや第二段階ROOTを再走査せず、既存の20期間histogram集計JSONから作る図。元Eのsource/config/helper、稼働中のwatcher、元ROOTを変更しない。描画はbatchと同じROOT5のstandalone Canvasを使う。

`ep8_data_001/ep8_cf_metrics.json` は7つの表示中央度と独立cent9 0–8、各plane −1/0–7の計144曲線を保存する。`plane=-1` はraw SE/MEを合算して一度規格化したinclusiveであり、8本のCFの平均ではない。0–10%だけの表示に加えて0–5%・5–10%も独立に残す。

- signalは背景差引き・純度補正前。
- H=[0.4,0.5) GeV/c、α=SE(H)/ME(H)、C=SE/(αME)。中央度はraw count合算後に規格化する。
- nativeは200bin、幅0.005 GeV/c。表示rebin4後の幅は0.02 GeV/c。
- 高k*窓は[0.5,0.6)、[0.6,0.8)、[0.8,1.0)、[0.6,1.0)、[0.5,1.0)。数値は窓内native countの比で、表示CF点の算術平均ではない。
- 上限1.0 GeV/cの外へ外挿しない。under/overflowは`countsIncludingFlow`へ別記する。
- 表示誤差はα固定のcount誤差。共通規格化、共有event/ME pairの相関を含む有意性・真の漸近値の判定は行わない。

0–60%の8planeはすべてHのSE/MEが正。SE(H)は14,423–15,257、ME(H)は2,750,687–2,959,186。C[0.8,1.0)は0.8990394–0.9234467、C[0.6,1.0)は0.9170951–0.9421395。inclusiveのC[0.6,1.0)=0.9302279。Hを基準とした観察結果であり、物理的plateauの認定ではない。

70–80%のEP0/2/3/4/5/7はHのSEがなくCF unavailable。EP1/6の高k*窓にSEがない場合の0.000は、規格化不能のNAと区別する。末梢binはH統計が少なく、α固定の誤差だけで強い比較をしない。exact countは`ep8_highk_windows.csv`にある。

処理済み保存eventは1,819,447,887、EP-validは1,736,875,306、EP-rejectedは82,572,581。全保存event数とEP有効なpair母集団を混同しない。

物理的な補正済み一次EP Ψ₁の周期2πを8分割する。sourceは`Float_t psi`、`Float_t half=psi/2`、丸めた幅`Float_t(pi/8)`を使い、halfの全体端点はstrictに除外する。図の45°刻みlabelはnominal。`canonicalPsi1.ranges`は実際のstored Float_tの各bin下端/上端をinclusiveで保存し、隣のFloat値が別binになることも確認する。Ψ₂の再構成を意味しない。

再現コマンド（既存結果を上書きしない新outdirが必要）：

```bash
python3 rootfile/lambda_ep8_cf_closure_report_20261120_001/prepare_ep8_cf.py \
  --settings rootfile/lambda_ep8_cf_closure_report_20261120_001/ep8_cf_settings.json \
  --output-dir rootfile/lambda_ep8_cf_closure_report_20261120_001/ep8_data_NEW \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_periodMix100.yaml
```

`--plot`を付けると`DrawEP8CF.C`でPDF/PNGを描く。figureも既存名を上書きしないので、再描画時は新settingsの`figureDirectory`へ新出力先を指定する。初回出力は`share/figure/20261120_lambda_ep8_cf_closure_report/`にある。元のprimary heatmapは保持し、読みやすさを整えた`ep8_highk_primary_heatmap_v2.pdf/png`を主資料に使う。v2の再現は新figure directoryを引数として次を実行する：

```bash
bash script/singularity_root_batch.sh \
  config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_periodMix100.yaml root4star -l -b -q \
  'rootfile/lambda_ep8_cf_closure_report_20261120_001/PolishEP8Heatmap.C+("TABLES_DIRECTORY","NEW_FIGURE_DIRECTORY")'
```

`ep8_data_001/ep8_completed.json`は入力SHA、20ROOT stat、111 protected production pinの前後不変と描画成功を記録する。元ROOTのpayload/tree entry読取は0。元の独立検証に対して0–60%の4窓×9曲線の36件を追加照合し、UF/OF込みの全cent9で8planeのraw合計がinclusiveと一致する。v2の入力TSVと元primary図の不変は`ep8_heatmap_polish_validation.json`に記録する。

この棚卸し・描画は、direct PicoDstとtree replayの同一prefix closure検証とは別の作業。その同値性の範囲や、mixing条件を変えた本番MEの解釈は親資料のclosure sectionを参照する。
