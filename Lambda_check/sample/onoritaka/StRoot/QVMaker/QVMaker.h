#ifndef _QVMAKER_H_
#define _QVMAKER_H_

#include "QV.h"
#include "ConstVar.h"
#include <iostream>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"
using namespace std;

class QVMaker {
	public:
		QVMaker();
		QVMaker( int i_mode );
		~QVMaker();

		void Init();
		void doCalibration( QV* qvec );
		void setProfiles( TProfile *hrecX, TProfile *hrecY, TProfile *hfltC, TProfile *hfltS, TH1F *hEP[3] );
		void LoadParameters( TFile* in );

	private:
		int m_CalibMode;

		// Recentering parameters
		float SMDXYmn         [NqvCent][NqvZvtx][NsubSMD][2];
		float SMDXYsg         [NqvCent][NqvZvtx][NsubSMD][2];
		float BBCXYmn[NordBBC][NqvCent][NqvZvtx][NsubBBC][2];
		float BBCXYsg[NordBBC][NqvCent][NqvZvtx][NsubBBC][2];
		float TPCXYmn[NordTPC][NqvCent][NqvZvtx][NsubTPC][2];
		float TPCXYsg[NordTPC][NqvCent][NqvZvtx][NsubTPC][2];
		float EMCXYmn[NordEMC][NqvCent][NqvZvtx][2];
		float EMCXYsg[NordEMC][NqvCent][NqvZvtx][2];
		float EPDXYmn[NordEPD][NqvCent][NqvZvtx][NsubEPD][2];
		float EPDXYsg[NordEPD][NqvCent][NqvZvtx][NsubEPD][2];
		float* get_RecPar(DETID fDet, int fSub, int fOrd, int fCent, int fVz, int fCS);

		// Flattening parameters
		static const int Nflt = 16;
		float SMDfltC         [NqvCent][NqvZvtx][NsubSMD][Nflt];
		float SMDfltS         [NqvCent][NqvZvtx][NsubSMD][Nflt];
		float BBCfltC[NordBBC][NqvCent][NqvZvtx][NsubBBC][Nflt];
		float BBCfltS[NordBBC][NqvCent][NqvZvtx][NsubBBC][Nflt];
		float TPCfltC[NordTPC][NqvCent][NqvZvtx][NsubTPC][Nflt];
		float TPCfltS[NordTPC][NqvCent][NqvZvtx][NsubTPC][Nflt];
		float EMCfltC[NordEMC][NqvCent][NqvZvtx][Nflt];
		float EMCfltS[NordEMC][NqvCent][NqvZvtx][Nflt];
		float EPDfltC[NordEPD][NqvCent][NqvZvtx][NsubEPD][Nflt];
		float EPDfltS[NordEPD][NqvCent][NqvZvtx][NsubEPD][Nflt];
		float* get_FltPar(DETID fDet, int fSub, int fOrd, int fCent, int fVz, int fCS);

		void doGainCorrection();

		void doNormalization( QV* qvec );
		void doRecentering  ( QV* qvec );
		void doFlattening   ( QV* qvec );

		TProfile *hrecX;
		TProfile *hrecY;
		TProfile *hfltC;
		TProfile *hfltS;
		TH1F *hEP[3];
};

#endif
