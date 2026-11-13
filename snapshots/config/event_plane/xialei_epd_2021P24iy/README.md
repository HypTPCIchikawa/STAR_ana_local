# Xialei EPD calibration: 2021 P24iy

このディレクトリは、Xialei さんの EPD 法で参照される recenter / shift 校正と、**校正の run index 順序を保った 388 run の表**を保存する。元ファイル、SHA256、run 表の照合、幾何コードの比較結果は [PROVENANCE.json](PROVENANCE.json) に記録した。

## 適用範囲

元のファイル取得スクリプトは `production_3p85GeV_fixedTarget_2021`、`P24iy`、`runnumber[]22158000-99999999` を明示している。元ディレクトリやファイル名の `3p0` という表記もこの production を指す。

この校正は production 全体を覆うものではない。保存 run 表は 388 run、最小 run は `22158032`、最大 run は `22179022`。**表に存在する run だけ**を使う。単に最小値と最大値の間に入っているだけでは適用できない。

既存の小 sample の run `22123040` は表に含まれず、この校正を適用できない。未知の run を別の run index に置き換えない。解析コードは coverage 不足、校正 SHA の不一致、必須 profile の欠落・bin 数不足、使用する校正 bin の空 entry / 非有限値を検出して停止する。全入力の coverage は実際の input list に対して別途確認する。

## Run index の根拠

次の run 表と `runs_2021.txt` は 388 run の全順序が一致することを確認した。

- Xialei の `tree_p24iy/StRoot_v2/run/run.h`
- jhan1 の `EP_P24iy/recenter/StRoot/run/run.h`
- 同 `shift/StRoot/run/run.h`
- 同 `shiftpar/StRoot/run/run.h`
- 上記 3 stage の `3p0GeV_local_P24iy_runmber_picodst.list`

校正生成コードは `run.h` の配列位置を run index として profile に `Fill` し、適用コードはその位置に対応する `runIndex+1` bin を読む。**run 表を sort し直さない。**

元の recenter profile 軸は centrality 9 bins `[-0.5, 8.5]` × run index 1000 bins `[-0.5, 999.5]`。shift profile はその間に shift term 20 bins `[-0.5, 19.5]` を持つ。

## 解析側の定義

EPD East の A/B はそれぞれ row 1–4 / 5–8、採用する角度は AB、すなわち row 1–8 の補正済み一次 EP `psi1shift_EPD[4]`。tile ごとに `RandomPointOnTile - primaryVertex` を使い、`nMIP >= 0.3`、weight を 2 で上限にし、weight 和で Q を規格化する。run と calibration centrality ごとの recenter 後、20 term の shift を適用する。

calibration centrality は全 primary track 数の閾値 `5, 10, 19, 33, 51, 77, 111, 155, 183` から計算する。CF の physics cent9 と別の量として保持する。原コードの event validity を引き継ぎ、両 TPC group と全 EPD group の一次・二次 Q の raw / recenter 条件も確認する。EP が無効な event は SE、ME、pool 保存のすべてから除外する。

mixing 用 8 bin は原 `Float_t` の角度変換と bin 幅の丸めを再現する。補正済み物理角度 `psi1` は独立に保存し、既存の TPC 二次 EP を上書きしない。1 bin の control でも同じ EP valid 条件を適用するため、8 bin と event 選択を揃えて比較できる。

## 幾何と乱数

Xialei の `StEpdGeom.cxx` と使用する STAR SL24y 版は、`ClassImp(StEpdGeom)` の有無以外に幾何・sampling 式の差がないことを確認した。row group の定義も一致する。

元の補正済み EP event cache は利用していない。設定 seed `20261112`、run ID、event ID から event ごとの seed を作る専用 `TRandom3` を使う。`StEpdGeom` は構築時の `gRandom` を保持するため、geometry の構築・sampling・破棄をこの専用 RNG の scope 内で行い、呼び出し元 RNG の状態を保つ。同一 event の再実行に対する再現性を持つが、元解析の event ごとの乱数列との一致を意味しない。

元 ROOT ファイルに書き込まない。校正を更新する際は対応する run 表、SHA256、生成 source の証拠も合わせて更新する。
