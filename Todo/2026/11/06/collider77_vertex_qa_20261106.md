# 7.7 GeV vertex QAの確認

状態：確認完了。7.7 GeV colliderでz=0に鋭いピークがないこと自体は異常ではなく、今回の選択後分布も幅広い衝突領域と整合する。これは頂点分布に関する限定した判断であり、P22ibのPID/purity等の全面的なQA合格を意味しない。稼働中のjob・解析設定・ライブラリは変更していない。

- [x] QA p.1左上のhVzは全14,046入力の選択前であり、bad-run/trigger/vertex/pileup/centrality選択後の図ではないことを確認。
- [x] primaryVertexのZをそのままfillし、FXT用のoffsetや表示座標変換がないことをコード監査で確認。
- [x] mode_vertex、Vz範囲±145 cm、Vr≤2 cm（中心0,0）、VPD cut明示無効を確認。vendored Run21 centrality表のVz範囲とも一致。
- [x] 7.7 GeV colliderの広いVz分布についてSTAR一次資料を確認。STAR Collaboration, Phys. Rev. C96,044904 (2017), Fig.2/本文p.4は7.7 GeVで中心付近が平坦になることを明記している。2010 BES-Iの例なので、2021 P22ibの個別QA合格を意味しない。
- [x] 選択後hCentralityVsVzのX projectionと元hVzを比較し、flow/積分・広がり・±55–60cm付近の局所構造を数値化した。
- [x] 比較図PDF/PNG・数値と説明を保存し、目視確認した。

閉じたpilot ROOTだけを読み取り、再event処理や解析jobの変更は行わない。hCentralityVsVzはevent/trigger/vertex/pileup・centrality計算後、cent-bin採否前にfillされる。今回の設定はcent9=0–8なので、projection積分をreconstructedEvents=3265と照合する。通常eventLedgerにVz/Vrがなく、eventごとの横方向分布はこのROOTから復元できない点を明記する。

参考：

- [STAR Collaboration論文（Fig.2とevent selection）](https://arxiv.org/pdf/1701.07065#page=4)
- [Horvat博士論文（Fig.4.2、本文p.39）](https://star.physics.yale.edu/Theses/Horvat_Thesis.pdf#page=54)：低energyでの広い縦方向分布、beam-pipe flangeに対応する±55cmのspikeと±75cm付近の構造を記載。今回のspikeの起源は直接測定ではなく、位置・選択前後比較に基づく推定として扱う。
- [準備Todo](collider77_preparation_20261106.md)
- [本番投入Todo](collider77_production_20261106.md)

追加図・数値の保存先：独立workspaceの`rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/vertex_qa_20261106_001/`。

## 実測結果

元hVzは14,046 Entries、可視範囲−150〜150 cmの積分10,171、underflow1,485、overflow2,390で閉じる。選択後projectionのEntries/積分は3,265、flow0でreconstructedEventsと一致し、全てのfine binでselected≤rawを確認した。細bin中心から計算した選択後の平均は−4.07 cm、RMS79.58 cmであり、この選択・可視範囲に条件づけた値である。colliderの幅広い分布で、0cmの狭いピークを要求する根拠はない。

選択前の−60〜−55cm付近は49.14 counts/cm（左右の隣接帯31.05/30.89）、+55〜+60cm付近は55.24 counts/cm（38.67/27.81）と局所増加がある。一方、選択後はそれぞれ12.19（12.95/11.33）、10.48（11.11/10.48）counts/cmで周囲と同程度となる。元図の鋭い側方spikeを正常な衝突領域のピークと解釈しない。既知のbeam-pipe構造の位置と整合するが、今回の入力で背景起源を直接同定したものではない。数値帯の実bin端・幅はsummary.jsonに記録した。

比較図はraw fine bin（0.75 cm）、選択後6cm bin、可視範囲の面積を各々1にしたshape比較、selected/rawのbinomial誤差付き割合を示す。割合は全event selectionの結果であり、vertex cut単独の効率ではない。6cmの端binは±145cm cutをまたぐため、bin中心±147cmの点は範囲外イベントを意味しない。

[比較図PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/vertex_qa_20261106_001/vertex_qa.pdf)／[比較図PNG](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/vertex_qa_20261106_001/vertex_qa.png)／[数値JSON](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/vertex_qa_20261106_001/summary.json)。ROOT5で閉じたpilot ROOTのhistogramだけを読み、以前のQA PDFは上書きせず追加作成した。
