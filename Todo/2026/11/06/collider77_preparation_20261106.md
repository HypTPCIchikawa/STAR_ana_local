# 7.7 GeV collider解析の準備

状態：ソフトウェア準備・少統計動作検証が完了。本番jobは未投入、P22ib用4He purity/校正は未確定。ユーザーが11月4日のTodoに沿った準備と本日のTodo作成を指示した。3.85 GeVのjobは終了済みだがhadd/検証/メール通知が稼働中のため、元のruntime/config/libraryを変更せず独立workspaceで進める。本番の全件job投入はこの準備作業には含めない。

対象：`sqrt(sNN)=7.7 GeV`、collider、2021年P22ib。解析はLambda–4He KFParticle。8次元`hNuclearPidPairJoint`を無効化し、候補cacheも作らないhist-only設定を準備する。通常event/候補台帳は標準出力として保持する。

## 実施項目

- [x] 現行hadd/通知の稼働状態と保護対象を確認する。
- [x] 必要な未追跡設定・ソースを含む独立workspaceを作る。library/build/ACLiC/tmp/job/log/err/rootfileを共有しない。
- [x] 2021 P22ibの正式なtrgsetupname、library、入力path/run/schema、catalog統計を確認する。
- [x] 7.7 GeV・2021年のcentrality parameter coverageと使用するmultiplicityを確認する。欠落時は別年やFXTへfallbackしない。
- [x] collider用Vz/Vr、trigger、KFのPV/幾何、rapidity、ME頂点範囲とcent bin、acceptance histogram範囲を調整する。
- [x] 7.7 GeV専用mainconf/analysis_infoと各YAMLを作成する。anaNameと出力先を独立させる。
- [x] 全核種をGit標準PID tableでpilot実行。R=2でのTOF切替、He4/d同一m²窓、R/rigidity vetoを明記。P24iy専用He4校正の適用ガードは保持し、無検証transferはしない。
- [ ] P22ibの十分な統計で4HeのPID係数/domain・purityを決める。今回He4候補1件では未確定。
- [x] 独立workspaceでbatch-matched buildを行い、元libraryへの参照が残らないことを確認する。
- [x] 少数の実入力を使って解析を実行し、event/Lambda/核種QA/ペア・通常histのfill、8D/候補cache不在を確認する。
- [x] 少統計のQA・出力量とcatalog統計から本番投入の設定・容量・未確認事項をまとめる。
- [x] 本番用job XML・手順README、再現設定と検証結果、日報を更新する。本番jobは未投入として区別する。

予定workspace：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib`。作成済み。15 libraryを新領域で再buildし、元libraryのRPATH/inode共有なしを確認。

## 初期catalog確認

bounded queryで`production_7p7GeV_2021`、P22ib、year2021、production library SL22b、ReversedFullFieldを確認した。最初の標本はrun22031042、`st_physics_22031042_raw_6500001.picoDst.root`、catalog7,824events/169,119,123bytes。全件入力には`filename~st_physics,filename!~adc`を明記する。productionのSL22bと、解析実行用libraryTagのABI/互換性は別途確認する。catalog値は実測PicoDst Entriesとは区別する。

[11月4日の調査・保留履歴](../04/collider77_preparation_20261104.md)を継承する。そこでのcatalog照会は未実行であり、データ不存在を意味しない。今回も検証済み/未確認を分けて記録する。Git commit/push、既存ROOT削除、実行中haddの変更は行わない。


## 実装・検証結果

- mainconf：`config/mainconf/main_auau7p7_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`（独立workspace内）。実行SL22b/ROOT5.34.38。analysis/centralityはrefmult、Run21 index41とbad-run/pileupを使用。
- 頂点±145 cm/Vr2 cm、MB OR810010/810020/810030/810040、VPD無効。mixing±145 cm/29 bins/depth20/cent9別。既存13.5 GeV設定もFXTであり、その頂点は使わない。
- 新providerのstrict trigger parser/gateを追加。MB各ID・OR・不正ID/空/重複/typo・reload、頂点境界/外側、旧Imp5互換性のintegration testがSL22bでPASS。変更sourceは独立領域のFemtoLambdaKfProvider.cxxだけ。
- Pilotの実PicoDst Entriesは4921+9125=14046。catalog申告7824+16887=24711とは異なるが、actual requested/input/completed14046と完全一致。2入力の20 branch構成も一致。全productionの実event総数は未測定。
- 14046入力、3265選択、Lambda1034、mixing clamp0。1786通常hist、通常ledger/counter一致、Sparse/候補cacheなし。He3候補27（signal SE12/ME154）、He4候補1（SE0/ME4、low-R/no-TOF）。38 QA histの内容・誤差が有限。
- manifest30511 unique URIsは不整合replicaの1 LFNを明示除外。RFF30493＋FullField18を保持。PilotはRFF2本でありFullFieldの実データ動作は未検証。
- 本番候補canonical XMLは`job/joblist/joblist_auau7p7_anaFemtoLambda_4He_KFParticle.xml`。input=filelist、max5、xrootd/local staging/ZIP。anaName/output stem/runmeta名を揃えた。生成済みだがsubmitは未実行。
- 出力量は1000 eventsで2,949,805 bytes、14046で3,711,750 bytes。2点modelのscenarioはROOT約40–70GB＋log/job最小割当約5–8GB、merge/stagingは追加。source統計や分割が未測定なので保証ではなく、投入直前に旧hadd領域を確保してquotaを再確認する。
- 元haddは01:45 JSTに10 checkpoints完了・継続を確認。元の通知/統合コード・mainconfなど6SHAは開始時と一致。元runtime/config/libraryへの変更なし。

[準備手順・全結果](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/preparation_20261106_001/README.md)／[QA図PDF（5ページ）](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/preparation_20261106_001/qa/pilot_14046_v2/collider77_pilot_qa.pdf)／[統合検証記録](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer_auau7p7_2021_P22ib/rootfile/auau7p7_anaFemtoLambda_4He_KFParticle/tests/preparation_20261106_001/validation.json)。

残作業はP22ibの4He PID/purity、より広いrun/fieldでのphysics QA、実際のSUMS分割と本番容量の確認。本番解析を検証済みとは扱わない。Git commit/pushなし、元ROOT削除なし。

## 続報：容量確認後の本番投入

上記のmax5/本番未投入は準備完了時点の記録として保持する。その後のユーザー指示により、08:12:46 JSTに隔離workspaceの本番投入をtmuxで開始した。analysis_infoのfilesPerJobを10にし、SUMS min6/max10・scratch20 GBを実simulationで検証した。物理設定・source/libraryは準備時と同じ。ID`36B7E875A19554DB711364BE0B072DA5`、最初の最大index5084で予定5085jobs。全件投入はまだ確認中。空き173.56 GBに対し、新解析85＋旧hadd60＋安全20＝165 GBを予約し、容量guardと完了後メールを設定した。以後の実投入・容量・状態は[本番Todo](collider77_production_20261106.md)に記録する。
