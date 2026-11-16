# Λ–4He：高k*のベースラインとME/k*計算の調査

状態：監査完了。低下の原因は未特定。復元org版p.18で0.8 GeV/c付近まで傾きがあるため、0.4–0.5の規格化は物理的plateauを確認した設定と扱わない。

- [x] org版p.18とL版のraw/補正後CFを0–1 GeV/cで確認
- [x] 提出時04BFのk*・四元運動量・混合条件をφ–d実装と照合
- [x] 不変量式・既知PRFの数値試験
- [x] Vz範囲外clamp、centrality合算、重み、paircap、両方向とeventIDを調査
- [x] 観測事実、除外できたバグ、未確定の原因と追加測定を区別して記録

既存PDF・ROOT・PID・production設定は変更せず、調査結果は新規ディレクトリcf_mixing_kstar_audit_20261102_001へ保存する。

## 結論

高k*の単調な低下は確認できた。一方、実投入04BFコードのk*計算、4Heの電荷補正、SE/ME共通処理、MEの基本的なfillに今回の低下を説明するバグは見つからなかった。**0.4–0.5 GeV/cは比較用の基準スケールであり、物理的plateauとして認定しない。原因を特定・検証するまでは最終CFの規格化根拠としない。**

- [図付き調査PDF（8頁、新規ファイル）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_mixing_kstar_audit_20261102_001/lambda_he4_mixing_kstar_audit_20261102.pdf)
- [数値・図・再現情報](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_mixing_kstar_audit_20261102_001/README.md)
- [検証](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_mixing_kstar_audit_20261102_001/validation.json)

### 確認できたこと

- 実提出04BF runtimeの主要6関数は現在と完全一致。全Maker差分は非Lambda用のgeneric処理のbad-run/pileup条件だけで、対象Lambda分岐はその前にreturnする。
- k*=|pA*−pB*|/2=|pA*|。ROOT5.34/38でLambda–4He・phi–d・Lambda–3He各20000ペア、計60000を既知PRF、不変量式、粒子交換、共通boostで検証してPASS。Lambda–4HeのPRF最大誤差7.55e−15 GeV/c、Float保存の差も最大1.07e−6 GeV/c。
- 4HeはpMomのrigidityを1回だけ2倍して質量3.72742でP4を作る。Lambdaは固定質量1.115683、再構成massは別変数。SE/MEとphi–dの共通StFemtoMakerで同じk*関数を使用する。phi–dのデータCFそのものは今回再解析していない。
- 実混合はbufferAll、両方向、同じcent9、Vz10bin（198–202cm）、pool/maxMixとも20。MEのeligible=attempted=filled=157569769、forward78787529、reverse78782240。cap、pair-cut、共有trackによる脱落0。
- eventIndexは各入力で増加し、現在イベントはmix後に保存する。共有track判定はeventIndexも比較する。mass classごとにmixを繰り返す重複はない。
- 懸念したVz範囲外clampは保存カウンタ0。粗いQAヒストの201–202.5cm binは境界202を跨ぐので、bin内全カウントを範囲外と数えてはいけない。このrunの原因としてclamp説は支持されない。

### データの切り分け

全て0.4–0.5を基準にした診断値。0.75の点は0.74–0.76の単一bin、0.7–0.8の値は区間のcount積分であり、同じ値ではない。

|中央度|0.75のraw|0.75のMethod5 L|cent9別ME露出を合わせる診断|
|---|---:|---:|---:|
|20–80%|0.88136±0.01276|0.87345±0.01279|0.87932±0.01273|
|0–20%|0.90792±0.00509|0.90570±0.00512|0.90640±0.00508|
|10–20%|0.88517±0.01034|0.88137±0.01038|0.88517±0.01034|

10–20%のraw0.7–0.8積分値は0.87625。細分cent9でも低下し、中央度の合算だけでは説明できない。cent9別ME exposure合わせもほとんど回復しない。この合わせ方は診断であり、物理補正には採用していない。導関数によるcount/規格化再利用の誤差は独立に点検済み。誤差は保存された統計情報の範囲で、event・shared-mixed-pair相関は含めない。

既存3He対照も調べた。0–20%、0.75では物理3He raw0.96128、同じ3Heを4He質量で計算0.95840、選別4He raw0.90766。対照にも傾きがあるが、質量仮説変更だけで4Heの全低下は再現しない。raw対照は純度補正しておらず、独立な無相関試料を証明するものではない。

L抽出はmass-histのwhole-bin窓、species対照はMakerの厳密なcandidate mass窓。両者の4He比の差は0.3–0.8で最大0.00102と、今回の低下より小さい。窓の違いは記録したままにした。

### 未解明の点と次の測定

今回の結果は、検証したk*式・基本ME処理にバグが見つからなかったという範囲であり、ME referenceが物理的に正しいと証明したものではない。理論のFSI成分が1に近いことだけで、実験SE/MEのイベント条件付き分布が等しくなるとは限らない。

候補/PIDを固定し、混合条件を1つずつ変えるcontrolled remixが次段階。cent9内部のより細かいmultiplicity・Vz、run/磁場、event-plane、保存則/flowと検出器・TOF受容の寄与を切り分ける必要がある。現在の投影TH1/TH2だけから任意のevent-level remixは作れず、Lambdaと核の両方の4元運動量・event classを含む候補記録が必要。今回は本番再処理を投入していない。規格化域を単に高kへ移すだけでは傾きの解決にならない。

既存のorg版・L-only v2レポートとROOTは保持。PID/Maker/config変更、production job、hadd、Git commit/pushは行っていない。

[最終照合・個人snapshot記録](mixing_audit_final_checks_20261102.json)。小さなコード・図・PDFは個人管理へ保存し、全hashを照合した。ROOT・大きな数値出力・ログ/バイナリは既存policyで除外し、解析側に保持する。
