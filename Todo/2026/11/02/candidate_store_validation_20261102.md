# Λ–4He候補保存・再混合の小統計検証（2026-11-02）

状態：実装・小統計検証・文書化を完了。

## 検証結果

| 項目 | 結果 |
|---|---|
| 入力 | 2021 P24iyのPicoDst 1ファイル、4,027イベント |
| acceptedイベントの保存・再生 | 3,858行、元の順序を維持 |
| 通常解析と保存追加解析 | ペア用TH1/2/3 1,675個＋Sparse 1個が完全一致 |
| 通常解析と候補だけの再生 | 同じ1,675個＋1個が完全一致 |
| 比較cell数 | 通常39,963,446（under/overflow含む）、Sparse 89,815 |
| 一致項目 | 軸・bin内容・誤差・Sumw2・Entries・統計量。許容差なし |
| 4He信号ペア | SE 9、ME 251（保存追加・元条件での再生とも一致） |
| Vz混合binを10→20へ変更 | SE 837個/19,981,638 cellは完全一致、4He ME251→235、native bin13個が変化 |
| 再生時 | PicoDst readerなし、LambdaProviderなし。既存StFemtoMakerの同じペア処理を使用 |
| ROOT5 synthetic試験 | 102 checks PASS（Float/P4/energy、空event、順序、不正payload、provenance、UUID、完了marker等） |
| 実Makerでの拒否試験 | PID変更／未許可のmixing変更をexit2で拒否。元cache不変 |

[検証図（1頁PDF）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_20261102_001/candidate_store_validation.pdf)。rebinは表示用のみ、CF規格化は行っていない。小統計で保存・再生方式の動作を確認したもので、高k*ベースラインの原因や最適な混合条件を決定したものではない。

## 保存内容・再生範囲

`femtoCandidateEvents` TTreeにacceptedイベントを順番に保存する。元のFemtoCandidateの20 Float、10 Int、UInt品質flagを直接保存し、energyも再計算しない。物理候補とpre-PID核種QA候補を分け、空collectionやvectorの順序も保持する。

イベントにはrun/event番号、元の入力番号、頂点、磁場、refMult/fxtMult、補正refMult/centWeight、cent9/16、track/TOF数、元のpsi2とvalid flagを保存する。再生時はイベント参照と娘track参照をstream内で一意な番号へ置き換える。

Λ・4He・3Heに加え、既存処理が集めるd/tも保存する。今回の物理候補はΛ62、4He250、3He887、d10,066、t1,803。物理候補なし726イベント、d/tのみ2,199イベントも保持した。d/tは直接ペアにしなくても元のpool登録条件に関わり、省略するとMEが変わる。QA候補211,297件、うち選択済み1,137件を保存した。

再生成するのは保存候補によるSE/MEのk*・mass・centrality/mT・PID pair QA等。PicoDst/helixに依存する再構成・track QA・SE merging QAは元の保存追加解析ROOTを参照する。再生ROOTのinputEventsはaccepted候補行数であり、元PicoDst入力4,027とは異なる。

## 混合条件・入力の扱い

元のPID・KF・核種選択・本番mainconfの設定値は変更していない。candidateStoreの指定がないmainconfでは機能は無効。selection provenance（参照先dEdx calibrationを含む）・species/schema・producer完了状態を読込前に確認する。mixing変更は専用設定で明示的に許可する。

候補入力は保存を行った元の各job ROOT、または順序付きROOT list。元TFileのUUIDを照合するため、候補treeをhadd/repackしたファイルは使わない。元条件のclosureにはjob境界でpoolをresetする。ファイルを跨ぐ再混合はcandidateStoreContinue mainconfで明示する。再生後のペアヒストグラムは通常どおり統合できる。

補助TPC Q2/psi2を別保存した。primary track、pT重み、候補track/娘を除外した未flattenのproxyであり、元のpairing psi2=-1は変更しない。valid3,857 / invalid1。EP bin>1の再生ではstored EPが必要で、invalidイベントは現状エラーとなる。較正済みEPや欠測処理は未検証。multiplicity/run/磁場等の情報も保存済みだが、それぞれの新しい混合bin方針の実装・物理比較は次段階。

candidate tree圧縮量8,934,587 bytes（非圧縮31,683,112 bytes）はtreeのみの値であり、出力ROOT全体のサイズや全統計への容量見積もりではない。

## 実装・再現

- [StFemtoMaker](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/StMaker/StFemtoMaker/StFemtoMaker.cxx)：既存ペア処理に保存・再生を接続。
- [中立な保存format](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/include/FemtoCandidateEventStore.h)／[設定・provenance](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/include/FemtoCandidateEventStoreConfig.h)。
- [保存mainconf](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWrite.yaml)／[元条件で再生](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreReplay.yaml)／[Vz20へ再混合](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreRemixVz.yaml)。
- [一括検証runner](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/tools/kfparticle_summary/validate_femto_candidate_store.py)／[保存・再生の説明](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/StMaker/StFemtoMaker/README.md)。
- [出力index・コマンド](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_20261102_001/README.md)／[主要検証記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_20261102_001/verification_003/validation.json)／[cache全行検査](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_20261102_001/cache_inspection.json)／[拒否試験](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_20261102_001/rejection_guard_tests.json)。

batch-matched SL24y / ROOT5.34/38でbuild。libStFemtoMaker.so SHA256: `b2dc1a6e1d869ca171f7880668c6e2ea761f5fad07692e2252b50f9053bda735`。元データrun22123040、入力SHA256: `3e300ad538179818845e74f571daf9af91986c0b55f3895b2f7af26925edd4b1`。

本検証はlocal実行。参照productionはSUMS `04BF2DAE47164B6F040E3A460E539289`（2021 P24iy）、mainconf `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml`、configlog `job/run/configlog/config_auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289.txt`。新しいSUMS jobは投入していない。

## チェックリスト

- [x] 既存ストアとFemtoCandidateの契約を確認。
- [x] raw FloatのP4/全候補情報・event情報・補助TPCQ2を保存。
- [x] 同じStFemtoMakerで再生（PicoDst/KF再構成なし）。
- [x] ROOT5 synthetic roundtrip・不正入力/provenance試験。
- [x] batch-matched build、実2021入力4,027イベントで基準/保存/再生を比較。
- [x] SE/MEのk*・mass・cent/mT・PID pair QA一致を確認。
- [x] 候補固定でVz binだけ変え、SE維持・ME変化を確認。
- [x] 複数streamで重複した元event番号を持つMEの照合。
- [x] 文書・日報・個人snapshotの最終保存。

全統計production投入・hadd・Git commit/pushは行っていない。

## 複数ファイルのevent参照の検証

実データから同じcent9/Vz binの2イベントを選び、1ファイルの入力番号[1,2]と、2ファイルの入力番号[1,1]を比較した。イベント/娘参照の意図した置換以外は20 Floatを含め元候補を固定。pool継続設定で再生し、5ヒストグラム・5,018 cellの軸・内容・誤差・Entries・統計量が完全一致した。4He信号ME=1、sampler eligible/attempted/filled=1/1/1、overlap veto=0が両方で一致。PicoDst/KF再構成は行っていない。

[複数ファイル検証](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_20261102_001/identity_fixture/replay_identity_comparison.json)／[全体の最終記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/candidate_store_20261102_001/completed.json)／[個人snapshot・最終確認](candidate_store_final_checks_20261102.json)。
