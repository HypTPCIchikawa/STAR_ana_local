#ifndef StBES2QaMaker_h
#define StBES2QaMaker_h

#include "StMaker.h"
#include "QVMaker/QV.h"
#include "QVMaker/QVMaker.h"
#include "StThreeVectorF.hh"
#include "SgPEvent/SgPEvent.h"
#include "SgPEvent/SgTrack.h"
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class TString;
class TH1F;
class TH2F;
class TH2D;
class TProfile;
class TProfile2D;
class StRefMultCorr;
//class StV0TofCorrection;
class TGraph;
class TF1;
class TVector3;

class TTree;
//class EventClass;

class StEpdEpFinder;

class StBES2QaMaker : public StMaker {
  public:
    StBES2QaMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StBES2QaMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    void    DeclareHistograms();
    void    WriteHistograms();

	void SetRunNumber(Int_t runnum){ mRunNumber = runnum; }
	void SetCalibrationMode( int mode );
		// 0= pass0: make gain correction factors of BBC-ADC 
		// 1= pass1: make recentering parameters
		// 2= pass2: make flattnining parameters
		// 3= analysis mode
	void SetDCAMax( float fDcaMax );
	void SetNHitsFitMax( int fNHitsFitMax );
	void SetZVertexMax( float fVzMax );
	void SetVertexRMax( float val ){ mVrAbsMax = val; }
	void SetVpdVzDif  ( float val ){ mVpdVzDif = val; }

	void SetTrgEffCorrection( int fTrgEff );
	void SetData( int flg );
		// 0= run10
		// 1= run11
		// 2= run14
		
	int GetRunNum(){ return mRunNumber; }
    
	void SetTopologicalCut( int fCut );
	void SetV0TofCorrection( int fCor );
	//void SetEvMixVtxShift( int flg );

	void SetRunRange(int min, int max){ 
		mRunmin  = min; 
		mRunmax  = max; 
		mRunbins = max - min;
		cout<<"StBES2QaMaker: SetRunRange("<<min<<","<<max<<")"<<endl; 
	}
	void SetFixedMode(int val){
		mFixedMode = val;
		cout<<"StBES2QaMaker: SetFixedMode("<<val<<")  0=collider 1=fixed-target"<<endl; 
	}
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
	int        mRunNumber;
	int        mZDCPEDRun;
	int        mCalibMode;
	int        mTrgSelection;
	float      mDcaMax;
	int        mNHitsFitMax;
	float      mVzAbsMax;
	float      mVrAbsMax;
	float      mVpdVzDif;
	int        mTrgEffCorrection;
	int        mTrkEffCorrection;
	float      mTrgEff;
	int        mData;  // 0=run10, 1=run11, 2=run14, 3=run16, 4=run18, 5=run19
	int        mTrkEP; // 0=primary, 1=global
	int        mEventCounter;
	int        mRunmin;
	int        mRunmax;
	int        mRunbins;
	int        mFixedMode;

	static const int ntrg = 6; // needs to be checked!

	StEpdEpFinder *mEpFinder;
	bool creatingWeights;

    
	Double_t refMult_pre;
	Float_t vz_pre;
	Float_t zdcAdcSumE_pre;
    
	int mTopologicalCut;
	int mEvMixVtxShift; // 0=off, 1=on

	float mDcaP [20];
	float mDcaPi[20];
	float mDcaV0[20];
	float mDcaDh[20];
	float mDL   [20];

	float mQvCh[4][4][2]; // [harmonic][sub x ch][x,y]
	float mQwCh[4][2];    // [sub x ch][mult,pTsum]

    TString    mOutName;
	StRefMultCorr *pRefMultCorr;
	QVMaker *pQVMaker;
	//StV0TofCorrection *pTofCorr;

	//tree
	//EventClass *pEve;
	TTree *ptree;

	static const int NCent = 9; // for vn measurements
	// 0-5%, 5-10%, 10-20%, 20-30%, 30-40%, 40-50%, 50-60%, 60-70%, 70-80%
	
	static Double_t mZDCSMDCenterex, mZDCSMDCenterey;           // ZDCSMD Beam Center
	static Double_t mZDCSMDCenterwx, mZDCSMDCenterwy;           // ZDCSMD Beam Center
	float mZDCSMDmap[2][2][8];

	float BbcTilePhi[2][18];
	Float_t GetBBCPhi(const Int_t eastWest, const Int_t tileId) const;
	void SetBBCTilesPhi();
	float GetBBCTilePhi(const Int_t eastwest, const Int_t tileId);
	void SetBBCAdc( float adcP[24], float adcT[18], int iep );
    
	Float_t ZDCSMD( StPicoEvent *pEv, int eastwest, int verthori, int strip ) const;
	Float_t ZDCSMD_GetPosition( int eastwest, int verthori, int strip );

	void CalibParamInit();
	Int_t LoadBBCalibParam( int runnum );
	Int_t LoadEPCalibParam( int runnum );

	void CalcQvSMD ( StPicoEvent *pEv, QV fQv[NsubSMD],          int fCent, int fVz );
	void CalcQvBBC ( StPicoEvent *pEv, QV fQv[NordBBC][NsubBBC], int fCent, int fVz );
	void CalcQvTPC ( StPicoDst   *pD,  QV fQv[NordTPC][NsubTPC], int fCent, int fVz );
	void CalcQvEPD ( StPicoDst   *pD,  QV fQv[NordEPD][NsubEPD], int fCent, int fVz );

	float SMDXYmn[NqvCent][NqvZvtx][NsubSMD][2];
	float SMDXYsg[NqvCent][NqvZvtx][NsubSMD][2];
	float BBCXYmn[NordBBC][NqvCent][NqvZvtx][NsubBBC][2];
	float BBCXYsg[NordBBC][NqvCent][NqvZvtx][NsubBBC][2];
	float TPCXYmn[NordTPC][NqvCent][NqvZvtx][NsubTPC][2];
	float TPCXYsg[NordTPC][NqvCent][NqvZvtx][NsubTPC][2];

	static const int Nflt = 16;
	float SMDfltC[NqvCent][NqvZvtx][NsubSMD][Nflt];
	float SMDfltS[NqvCent][NqvZvtx][NsubSMD][Nflt];
	float BBCfltC[NordBBC][NqvCent][NqvZvtx][NsubBBC][Nflt];
	float BBCfltS[NordBBC][NqvCent][NqvZvtx][NsubBBC][Nflt];
	float TPCfltC[NordTPC][NqvCent][NqvZvtx][NsubTPC][Nflt];
	float TPCfltS[NordTPC][NqvCent][NqvZvtx][NsubTPC][Nflt];

	bool IsTofMatch( double& refMult, unsigned short& tofMult, int runnum );
	//float getTofBeta(StPicoTrack const* const trk, StThreeVectorF const& vtx, float BField) const;
	float getTofBeta(StPicoTrack const* const trk, TVector3 const& vtx, float BField) const;

	//int GetVnCentralityBin( float refMult, int fVz, int ftrg );
	//int GetCentralityBin  ( float refMult, int fVz, int ftrg );
	//int GetZvtxBin        ( float zvtx    );

	bool isGoodPTrack( StPicoTrack *trk );
	bool isGoodGTrack( StPicoTrack *trk );
	bool isGoodTrackEP( StPicoTrack *trk );

	enum PAIRTYPE { REAL, MIX };
	void MakeV0Pair(SgPEvent& eveA, SgPEvent& eveB, int fCen, int fCH, PAIRTYPE fPair );

	bool isProtonTof(float fm2);
	bool isPionTof  (float fm2, float fp);
	//bool TofCorrM2Cut(StThreeVectorF& pvtx, StThreeVectorF& v0vtx, StLorentzVectorD& v0p, SgTrack& trk, StPhysicalHelixD& helix, int ipart );

	SgPEvent pEvt [2];
	SgPEvent piEvt[2];
	double mAch;
	double mAk;

	void TrackLoop(int& myrefMult, int& mygrefMult);


	// event info.
	// ___________________
	// vertex
	TH1F *hVz;
	TH1F *hdVz;
	TH2F *hVxVy;
	TH2F *hVxVyCut;
	TH2F *hVzVpdz;
	TH2F *hVzVr;

	TH1F *hVz_VrCut;

	// multiplicity
	TH1F *hrefMult;
	TH1F *hgrefMult;
	TH1F *htofMult;
	TH2F *hrefMultvsTofMult;
	TH2F *hrefMultvsTofMatch;
	TH2F *hrefMultvsETofMult;

	TH1F *hTrg;
	TH1F *hVtxRank[2];
	TH1F *hnTofT0;

	TH2F *hBBCEvsW;
	TH2F *hZDCEvsW;
	TH2F *hBBCvsMult [2];
	TH2F *hZDCvsMult [2];
	TH2F *hBBCvsGMult[2];
	TH2F *hZDCvsGMult[2];
	TH2F *hBBCvsTMult[2];
	TH2F *hBBCAdc2D;
	TH2F *hZDCAdc2D;
	TH2F *hBBCAdc2DHi;
	TH2F *hrefMultvsBBCRingE[4];
	TH2F *hrefMultvsBBCRingW[4];

	TH1F *hZdcRate[3];
	TH1F *hBbcRate[3];
	TH1F *hBGRate;

	TH1F *BBCgc;
	TH2D *hZDCSMDped;
	TH2F *hZDCAdc   [2];
	TH2F *hZDCAdcCor[2];

	TH2F *hVzvsrefMult;
	TH2F *hVzvsEpdMipE;
	TH2F *hVzvsEpdMipW;
	TProfile *hprVzvsrefMult;
	TProfile *hprVzvsEpdMipE;
	TProfile *hprVzvsEpdMipW;

	TH2F *hMtdHit;
	TH2F *hBTowHit;

	// track info.
	// QA mode
	TH1F *hNhitsFit;
	TH1F *hNhitsMax;
	TH1F *hgDca;
	TH1F *hDca;
	TH1F *hDcaXY;
	TH1F *hDcaZ;
	TH1F *hPPt[2];
	TH1F *hGPt[2];
	TH1F *hPdEdx;
	TH1F *hGdEdx;
	TH2F *hPPhivsEta[2];
	TH2F *hGPhivsEta[2];
	TH2F *hGPhivsEta2[2];
	TH2F *hPPhivsEtaTof[2];
	TH2F *hPPhivsEtaMtd;
	TH2F *hPPhivsEtaBemc;
	TH2F *hPdEdxvsP;
	TH2F *hGdEdxvsP;
	TH2F *hGdEdxvsP_PosPion3orii;
	TH2F *hGdEdxvsP_PosPion2orii;
	TH2F *hGdEdxvsP_NegPion3orii;
	TH2F *hGdEdxvsP_NegPion2orii;
	TH2F *hGdEdxvsP_Prton3orii;
	TH2F *hGdEdxvsP_Prton2orii;
	TH2F *hGdEdxvsProton;
	TH2F *hGdEdxvsPion;
	TH1F *hRatioPrGlTrk;

	// pid check
	TH2F *hM2vsP;
	TH2F *hM2vsPdca;
	TH2F *hDcavsP[6];
	TH2F *hNsigmaP[4];
	TH2F *hBetaIvsP;
	TH2F *hBetaIvsPe;
	TH2F *hETOFxy;

	// Lambda
	TH1F *hLambdaMass;
	//TH1F *hLambdaMass[16];
	//TH1F *hLambdaMassTopoCut[0];
	TH1F *hLambdaMassTopoCut;
	TH1F *hLambdaMassTopoCutpos;
	TH1F *hLambdaMassTopoCutneg;
	TH2F *hLambdaMassdphiAB;
	TH2F *hLambdaMassdphiAL;

	TH1F *hLambdaMassTopoCutCent[10];
	TH1F *hLambdaMassTopoCutCentpos[10];
	TH1F *hLambdaMassTopoCutCentneg[10];

	TH1F *hLambdaMassTopoCutPt[5];
	TH1F *hLambdaMassTopoCutPtpos[5];
	TH1F *hLambdaMassTopoCutPtneg[5];

	TH1F *hLambdaMassTopoCutRap[5];
	TH1F *hLambdaMassTopoCutRappos[5];
	TH1F *hLambdaMassTopoCutRapneg[5];
	TH1F *hLambdaMassTopoCutCPR[4][4][4];
	//TH1F *hLambdaMassTopoCutPtDtail[8];
	//TH1F *hLambdaMassPtCutEtaDtail[6];
	//TH1F *hLambdaMassTopoCutiCent[8];
	TH2D *hLambdaMassCPR[10][4][4];
        TH2F *hPhivsEtaLambda;
        TH2F *hPtvsEtaLambda;
        
    //EPD(orii)
    TH2D *qvHist[NordEPD][24][8];
    TProfile *qvProf[NordEPD][24][2];
    TProfile *qvFProf[NordEPD][24][10];
	TH2F *hPsivsPsi[2][24][24];
	TH2F *hLambdarapvspt;
	TProfile *hResCos[2][24][24];
	TProfile *hResSin[2][24][24];
	TProfile *hResCos2[24][24];
	TH1F *hLambdaphi[9];
	TProfile *hv1vsMinv[9];
	TProfile *hv1vsMinvpos;
	TProfile *hv1vsMinvneg;

	TProfile *hv1sinvsMinv[9];
	TProfile *hv1sinvsMinvpos;
	TProfile *hv1sinvsMinvneg;

	TProfile *hv1vsMinvCent[10];
	TProfile *hv1vsMinvCentpos[10];
	TProfile *hv1vsMinvCentneg[10];

	TProfile *hv1sinvsMinvCent[10];
	TProfile *hv1sinvsMinvCentpos[10];
	TProfile *hv1sinvsMinvCentneg[10];

	TProfile *hv1vsMinvPt[5];
	TProfile *hv1vsMinvPtpos[5];
	TProfile *hv1vsMinvPtneg[5];

	TProfile *hv1sinvsMinvPt[5];
	TProfile *hv1sinvsMinvPtpos[5];
	TProfile *hv1sinvsMinvPtneg[5];

	TProfile *hv1vsMinvRap[5];
	TProfile *hv1vsMinvRappos[5];
	TProfile *hv1vsMinvRapneg[5];
	TProfile *hv1vsMinvCPR[4][4][4];

	TProfile *hv1sinvsMinvRap[5];
	TProfile *hv1sinvsMinvRappos[5];
	TProfile *hv1sinvsMinvRapneg[5];

	TProfile *hLambdaphivspt[9];
	TProfile *hLambdaphivseta[9];
	TProfile *hLambdaphivscen[9];

	TH1F *hPrPiphi[9];
	TH1F *hLamdphi[9];
	TH1F *hLamNdphi[9];
	TH1F *hProLamdphi[9];
	TH1F *ppqqdphi;
	TH1F *hLamCosdphi[9];
	TH1F *hLamSindphi[9];
	TProfile2D *v1dphi;
	TProfile2D *v1sindphi;
	TProfile2D *v1dphiAL;
	TProfile2D *v1sindphiAL;
	TProfile2D *GlobalPoldphi;
	TProfile2D *GlobalPolcosdphi;
	TProfile2D *GlobalPoldphiAL;
	TProfile2D *GlobalPolcosdphiAL;
	TProfile2D *ProductionPoldphi;
	TProfile2D *ProductionPolcosdphi;
	TProfile2D *ProductionPoldphiAL;
	TProfile2D *ProductionPolcosdphiAL;
	TProfile2D *ProductionPolPt;
	TH2F *hAaLamvsAB;
	TProfile *hpolsinvsMinv[9];
	TProfile *hpolsinvsMinvpos;
	TProfile *hpolsinvsMinvneg;
	TProfile *hpolcosvsMinv[9];
	TProfile *hpolcosvsMinvpos;
	TProfile *hpolcosvsMinvneg;
	TProfile *hpolcosvsMinvCent[10];
	TProfile *hpolcosvsMinvCentpos[10];
	TProfile *hpolcosvsMinvCentneg[10];
	TProfile *hpolcosvsMinvPt[5];
	TProfile *hpolcosvsMinvPtpos[5];
	TProfile *hpolcosvsMinvPtneg[5];
	TProfile *hpolcosvsMinvRap[5];
	TProfile *hpolcosvsMinvRappos[5];
	TProfile *hpolcosvsMinvRapneg[5];
//	TProfile *hpolsinvsMinvCent[10];
//	TProfile *hpolsinvsMinvCentpos[10];
//	TProfile *hpolsinvsMinvCentneg[10];
	TProfile *hpolsinvsMinvPt[5];
	TProfile *hpolsinvsMinvPtpos[5];
	TProfile *hpolsinvsMinvPtneg[5];
	TProfile *hpolsinvsMinvRap[5];
	TProfile *hpolsinvsMinvRappos[5];
	TProfile *hpolsinvsMinvRapneg[5];
	TProfile *hplpolsinvsMinv[9];
	TProfile *hplpolsinvsMinvpos;
	TProfile *hplpolsinvsMinvneg;
	TProfile *hplpolcosvsMinv[9];
	TProfile *hplpolcosvsMinvpos;
	TProfile *hplpolcosvsMinvneg;
	TProfile *hplpolcosvsMinvCent[10];
	TProfile *hplpolcosvsMinvCentpos[10];
	TProfile *hplpolcosvsMinvCentneg[10];
	TProfile *hplpolcosvsMinvPt[5];
	TProfile *hplpolcosvsMinvPtpos[5];
	TProfile *hplpolcosvsMinvPtneg[5];
	TProfile *hplpolcosvsMinvRap[5];
	TProfile *hplpolcosvsMinvRappos[5];
	TProfile *hplpolcosvsMinvRapneg[5];
//	TProfile *hplpolsinvsMinvCent[10];
//	TProfile *hplpolsinvsMinvCentpos[10];
//	TProfile *hplpolsinvsMinvCentneg[10];
	TProfile *hplpolsinvsMinvPt[5];
	TProfile *hplpolsinvsMinvPtpos[5];
	TProfile *hplpolsinvsMinvPtneg[5];
	TProfile *hplpolsinvsMinvRap[5];
	TProfile *hplpolsinvsMinvRappos[5];
	TProfile *hplpolsinvsMinvRapneg[5];

	TProfile2D *v1dphiALCPR[10][4][4];
    TProfile2D *GlobalPoldphiALCPR[10][4][4];
	TProfile2D *ProductionPoldphiLRCPR[10][4][4]; 
	
	TProfile *hLambdaphivsptCent[9][3];
	TProfile *hLambdaphivsetaCent[9][3];
	TH2F *hPsivsPsiCent[2][24][24][3];
	TProfile *hChvspt[9][2]; 
	TProfile *hChvseta[9][2]; 
	TProfile *hChvsCent[9][2]; 
	TProfile *hChvsptCent[9][2][3];
	TProfile *hChvsetaCent[9][2][3];
    float Prm[NordEPD][24][2][10][2];
    float Prmf[NordEPD][24][10][16];
    float qQpsi[2][24];
	int imi;

	// EP 
	TH1F *hBBCEP[NordBBC][NsubBBC][NqvCent][3];
	TH1F *hSMDEP[NsubSMD][NqvCent][3];
	TH1F *hTPCEP[NordTPC][NsubTPC][NqvCent][3];
	TH1F *hEPDEP[NordEPD][NsubEPD][NqvCent][3];
	
	TH2F *hQSMD2D[NsubSMD][2];
	TH1F *hQBBC[2][NsubBBC][NCent];
	TH1F *hQTPC[3][NsubTPC][NCent];

	// EP calibration
	TProfile *hBbcAdc;
	TProfile *hBbcAdcCor;
	TProfile *hSMDrecX[NsubSMD];
	TProfile *hSMDrecY[NsubSMD];
	TProfile *hSMDfltC[NsubSMD];
	TProfile *hSMDfltS[NsubSMD];
	TProfile *hBBCrecX[NordBBC][NsubBBC];
	TProfile *hBBCrecY[NordBBC][NsubBBC];
	TProfile *hBBCfltC[NordBBC][NsubBBC];
	TProfile *hBBCfltS[NordBBC][NsubBBC];
	TProfile *hTPCrecX[NordTPC][NsubTPC];
	TProfile *hTPCrecY[NordTPC][NsubTPC];
	TProfile *hTPCfltC[NordTPC][NsubTPC];
	TProfile *hTPCfltS[NordTPC][NsubTPC];
	TProfile *hEPDrecX[NordEPD][NsubEPD];
	TProfile *hEPDrecY[NordEPD][NsubEPD];
	TProfile *hEPDfltC[NordEPD][NsubEPD];
	TProfile *hEPDfltS[NordEPD][NsubEPD];

	// EP correlation
	TProfile *hSMDEPCrrvsMult[NordSMD][2];
	TProfile *hBBCEPCrrvsMult[NordBBC][4];
                    
	//TProfile *APTvsPtBin[2][NCent];

	// only for FXT31
	TProfile *hQvCor[2][4];
	TProfile *hRunvsQvPar[2][2][10][2];
	TProfile *hInRunvsQvPar[2][2][10][2];

	TH1F *hNEpdHits;
	TH1F *hEpdHitAdc;
	TH2F *hEpdEtaPhiEast;
	TH2F *hEpdEtaPhiWest;
	TH2F *hrefMultvsEpdRingE[16];
	TH2F *hrefMultvsEpdRingW[16];
	TH2F *hrefMultvsEpdRingEraw[16];
	TH2F *hrefMultvsEpdRingWraw[16];
	TH2F *hEpdEvsW;

	// run dependece
	TProfile *hRunidvsrefMult;
	TProfile *hRunidvsgrefMult;
	TProfile *hRunidvstofMult;
	TProfile *hRunidvstofMatch;
	TProfile *hRunidvsBemcMatch;
	TProfile *hRunidvsEpdMipE;
	TProfile *hRunidvsEpdMipW;
	TProfile *hRunidvsEtofMult;
	TProfile *hRunidvsVx;
	TProfile *hRunidvsVy;
	TProfile *hRunidvsVz;
	TProfile *hRunidvsRank;
	TProfile *hRunidvsZdcx;
	TProfile *hRunidvsBbcx;
	TProfile *hRunidvsBgrate;
	TProfile *hRunidvsTpcQx[2];
	TProfile *hRunidvsTpcQy[2];
	TProfile *hRunidvsEpdQx[2][2];
	TProfile *hRunidvsEpdQy[2][2];
	TH1F     *hRunidvsAllEvt;
	TH1F     *hRunidvsGoodEvt[4];

	TProfile *hRunidvsPt;
	TProfile *hRunidvsEta;
	TProfile *hRunidvsPhi;
	TProfile *hRunidvsNchp; 
	TProfile *hRunidvsNchm;
	TProfile *hRunidvsNproton;
	TProfile *hRunidvsNprotonbar;
	TProfile *hRunidvsNhits;
	TProfile *hRunidvsDedx;
	TProfile *hRunidvsDca;
	TH2F     *hRunidvsEtaDist;
	TH1F *hPt;
	TH1F *hEta;
	TH1F *hPhi;
	TH1F *hNhits;
	TH1F *hDedx;
	TH1F *hDcaTmp;

	TH1F *hVzAll;

	TProfile2D *hAveDedx2D[2];
	TProfile *hAveDedxPhi[2];
	TH2F *hDedxPhi[2];

    ClassDef(StBES2QaMaker, 1)
};

#endif
