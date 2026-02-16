#include "QVMaker.h"

// constructor
QVMaker::QVMaker() : m_CalibMode(3) {
	Init();
}

QVMaker::QVMaker( int i_mode ) : m_CalibMode(i_mode) {
	Init();
}

//destructor
QVMaker::~QVMaker(){
}

void QVMaker::Init(){

	fill_n( &(SMDXYmn[0][0][0][0]),    sizeof(SMDXYmn)/sizeof(float), 0 );
	fill_n( &(SMDXYsg[0][0][0][0]),    sizeof(SMDXYsg)/sizeof(float), 1 );
	fill_n( &(BBCXYmn[0][0][0][0][0]), sizeof(BBCXYmn)/sizeof(float), 0 );
	fill_n( &(BBCXYsg[0][0][0][0][0]), sizeof(BBCXYsg)/sizeof(float), 1 );
	fill_n( &(TPCXYmn[0][0][0][0][0]), sizeof(TPCXYmn)/sizeof(float), 0 );
	fill_n( &(TPCXYsg[0][0][0][0][0]), sizeof(TPCXYsg)/sizeof(float), 1 );
	fill_n( &(EMCXYmn[0][0][0][0]),    sizeof(EMCXYmn)/sizeof(float), 0 );
	fill_n( &(EMCXYsg[0][0][0][0]),    sizeof(EMCXYsg)/sizeof(float), 1 );
	fill_n( &(EPDXYmn[0][0][0][0][0]), sizeof(EPDXYmn)/sizeof(float), 0 );
	fill_n( &(EPDXYsg[0][0][0][0][0]), sizeof(EPDXYsg)/sizeof(float), 1 );

	fill_n( &(SMDfltC[0][0][0][0]),    sizeof(SMDfltC)/sizeof(float), 0 );
	fill_n( &(SMDfltS[0][0][0][0]),    sizeof(SMDfltS)/sizeof(float), 0 );
	fill_n( &(BBCfltC[0][0][0][0][0]), sizeof(BBCfltC)/sizeof(float), 0 );
	fill_n( &(BBCfltS[0][0][0][0][0]), sizeof(BBCfltS)/sizeof(float), 0 );
	fill_n( &(TPCfltC[0][0][0][0][0]), sizeof(TPCfltC)/sizeof(float), 0 );
	fill_n( &(TPCfltS[0][0][0][0][0]), sizeof(TPCfltS)/sizeof(float), 0 );
	fill_n( &(EMCfltC[0][0][0][0]),    sizeof(EMCfltC)/sizeof(float), 0 );
	fill_n( &(EMCfltS[0][0][0][0]),    sizeof(EMCfltS)/sizeof(float), 0 );
	fill_n( &(EPDfltC[0][0][0][0][0]), sizeof(EPDfltC)/sizeof(float), 0 );
	fill_n( &(EPDfltS[0][0][0][0][0]), sizeof(EPDfltS)/sizeof(float), 0 );

	hrecX  = NULL;
	hrecY  = NULL;
	hfltC  = NULL;
	hfltS  = NULL;
	hEP[0] = NULL;
	hEP[1] = NULL;
	hEP[2] = NULL;
}

void QVMaker::doCalibration( QV* qvec ){

	doNormalization( qvec );
	doRecentering  ( qvec );
	doFlattening   ( qvec );
}

void QVMaker::doNormalization( QV* qvec ){

	if( qvec->get_Qw()>0 ){

		qvec->set_Qx( qvec->get_Qx()/sqrt(qvec->get_Qw()) );
		qvec->set_Qy( qvec->get_Qy()/sqrt(qvec->get_Qw()) );
		qvec->set_Psi( atan2( qvec->get_Qy(), qvec->get_Qx() ) / (qvec->get_Ord()+1.0) );

		if( m_CalibMode>0 && m_CalibMode<3 ){
			int id = 1 + qvec->get_Zvtx() + qvec->get_Cent()*NqvZvtx;
			hrecX->Fill( id, qvec->get_Qx() );
			hrecY->Fill( id, qvec->get_Qy() );
		}
		if( m_CalibMode==3 ) hEP[0]->Fill( qvec->get_Psi() );
	}else{
		qvec->set_Qx( -9999.0 );
		qvec->set_Qy( -9999.0 );
		qvec->set_Qabs( -9999.0 );
	}
}

void QVMaker::doRecentering( QV* qvec ){

	float* recMn = get_RecPar( qvec->get_Det(), qvec->get_Sub(), qvec->get_Ord(), qvec->get_Cent(), qvec->get_Zvtx(), 0); 
	float* recSg = get_RecPar( qvec->get_Det(), qvec->get_Sub(), qvec->get_Ord(), qvec->get_Cent(), qvec->get_Zvtx(), 1); 

	if( qvec->get_Qw()>0 && recSg[0]!=0 && recSg[1]!=0 ){
		qvec->set_Qx( (qvec->get_Qx()-recMn[0])/recSg[0] );
		qvec->set_Qy( (qvec->get_Qy()-recMn[1])/recSg[1] );
		qvec->set_Psi( atan2( qvec->get_Qy(), qvec->get_Qx() ) / (qvec->get_Ord()+1.0) );
		qvec->set_Qabs( sqrt(qvec->get_Qx()*qvec->get_Qx()+qvec->get_Qy()*qvec->get_Qy()) );
		if( m_CalibMode==3 ) hEP[1]->Fill( qvec->get_Psi() );
	}else{
		qvec->set_Qw( -9999.0 );
	}
}

void QVMaker::doFlattening( QV* qvec ){

	float* fltC = get_FltPar( qvec->get_Det(), qvec->get_Sub(), qvec->get_Ord(), qvec->get_Cent(), qvec->get_Zvtx(), 0 );
	float* fltS = get_FltPar( qvec->get_Det(), qvec->get_Sub(), qvec->get_Ord(), qvec->get_Cent(), qvec->get_Zvtx(), 1 );

	if( qvec->get_Qw()>0 ){

		float dPsi = 0.0;
		float nPsi = ( qvec->get_Ord() + 1.0 ) * qvec->get_Psi();
		for( int ik=0; ik<Nflt; ik++ ){

			float cterm = cos( (ik+1.0)*nPsi );
			float sterm = sin( (ik+1.0)*nPsi );
			if( m_CalibMode==2 ){
				int id = 1 + ik + qvec->get_Zvtx()*Nflt + qvec->get_Cent()*Nflt*NqvZvtx;
				hfltC->Fill( id, cterm );
				hfltS->Fill( id, sterm );
			}

			float fCos =  2.0 * fltC[ik] / (ik+1.0);
			float fSin = -2.0 * fltS[ik] / (ik+1.0);
			dPsi += fSin*cterm + fCos*sterm;
		}
		nPsi += dPsi;
		qvec->set_Psi( atan2( sin(nPsi), cos(nPsi) ) / (qvec->get_Ord()+1.0) );
		if( m_CalibMode==3 ) hEP[2]->Fill( qvec->get_Psi() );
	}else{
		qvec->set_Psi( -9999. );
	}

}

float* QVMaker::get_FltPar(DETID fDet, int fSub, int fOrd, int fCent, int fVz, int fCS){

	float* ptr;
	if( fDet==dZDC ){ //ZDC-SMD
		if( fCS==0 ) ptr = SMDfltC[fCent][fVz][fSub];
		else         ptr = SMDfltS[fCent][fVz][fSub];

	}else if( fDet==dBBC ){ //BBC
		if( fCS==0 ) ptr = BBCfltC[fOrd][fCent][fVz][fSub];
		else         ptr = BBCfltS[fOrd][fCent][fVz][fSub];

	}else if( fDet==dTPC ){ //TPC
		if( fCS==0 ) ptr = TPCfltC[fOrd][fCent][fVz][fSub];
		else         ptr = TPCfltS[fOrd][fCent][fVz][fSub];

	}else if( fDet==dEPD ){ //EPD
		if( fCS==0 ) ptr = EPDfltC[fOrd][fCent][fVz][fSub];
		else         ptr = EPDfltS[fOrd][fCent][fVz][fSub];


	}else{ //EEMC
		if( fCS==0 ) ptr = EMCfltC[fOrd][fCent][fVz];
		else         ptr = EMCfltS[fOrd][fCent][fVz];
	}

	return ptr;
}

float* QVMaker::get_RecPar(DETID fDet, int fSub, int fOrd, int fCent, int fVz, int fCS){

	float* ptr;
	if( fDet==dZDC ){ //ZDC-SMD
		if( fCS==0 ) ptr = SMDXYmn[fCent][fVz][fSub];
		else         ptr = SMDXYsg[fCent][fVz][fSub];

	}else if( fDet==dBBC ){ //BBC
		if( fCS==0 ) ptr = BBCXYmn[fOrd][fCent][fVz][fSub];
		else         ptr = BBCXYsg[fOrd][fCent][fVz][fSub];

	}else if( fDet==dTPC ){ //TPC
		if( fCS==0 ) ptr = TPCXYmn[fOrd][fCent][fVz][fSub];
		else         ptr = TPCXYsg[fOrd][fCent][fVz][fSub];

	}else if( fDet==dEPD ){ //EPD
		if( fCS==0 ) ptr = EPDXYmn[fOrd][fCent][fVz][fSub];
		else         ptr = EPDXYsg[fOrd][fCent][fVz][fSub];

	}else{ //EEMC
		if( fCS==0 ) ptr = EMCXYmn[fOrd][fCent][fVz];
		else         ptr = EMCXYsg[fOrd][fCent][fVz];
	}

	return ptr;
}

void QVMaker::LoadParameters( TFile* in ){

	if( !in ){ cout << endl << "Could not open a calibration file !!!!!" << endl << endl; return; }

	// SMD
	// ------------------------------------------------------------------------
	// Recentering parameters
	TProfile *SMDX[NsubSMD];
	TProfile *SMDY[NsubSMD];
	for( int iep=0; iep<NsubSMD; iep++ ){
   		SMDX[iep] = (TProfile*)in->Get( Form("hSMDrecX_%d",iep) );
   		SMDY[iep] = (TProfile*)in->Get( Form("hSMDrecY_%d",iep) );
		if( SMDX[iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get SMDX!"<< endl; return; }
		if( SMDY[iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get SMDY!"<< endl; return; }
	}
	for( int ict=0; ict<NqvCent; ict++ ){
		for( int iz=0; iz<NqvZvtx; iz++ ){
			int id = 1+iz+ict*NqvZvtx;
			for( int iep=0; iep<NsubSMD; iep++ ){
				SMDXYmn[ict][iz][iep][0] = SMDX[iep]->GetBinContent(id);
				SMDXYsg[ict][iz][iep][0] = SMDX[iep]->GetBinError  (id);
				SMDXYmn[ict][iz][iep][1] = SMDY[iep]->GetBinContent(id);
				SMDXYsg[ict][iz][iep][1] = SMDY[iep]->GetBinError  (id);
				//cout << SMDXYmn[ict][iz][iep][0] <<" "<< SMDXYsg[ict][iz][iep][0] <<" "<<  SMDXYmn[ict][iz][iep][1] <<" "<< SMDXYsg[ict][iz][iep][1] << endl;
	}	}	}

	// Flattening parameters
	if( m_CalibMode==3 ){
		TProfile *fltSMDC[NsubSMD];
		TProfile *fltSMDS[NsubSMD];
		for( int iep=0; iep<NsubSMD; iep++ ){
   			fltSMDC[iep] = (TProfile*)in->Get( Form("hSMDfltC_%d",iep) );
  			fltSMDS[iep] = (TProfile*)in->Get( Form("hSMDfltS_%d",iep) );
			if( fltSMDC[iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltSMDC!"<< endl; return; }
			if( fltSMDS[iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltSMDS!"<< endl; return; }
		}
		for( int iep=0; iep<NsubSMD; iep++ ){
			for( int ict=0; ict<NqvCent; ict++ ){
				for( int iz=0; iz<NqvZvtx; iz++ ){
					for( int ik=0; ik<Nflt; ik++ ){
						int id = 1 + ik+iz*Nflt+ict*Nflt*NqvZvtx;
						SMDfltC[ict][iz][iep][ik] = fltSMDC[iep]->GetBinContent(id);
						SMDfltS[ict][iz][iep][ik] = fltSMDS[iep]->GetBinContent(id);
						//cout <<" "<<ict<<" "<<iz<<" "<< SMDfltC[ict][iz][iep][ik] <<" "<< SMDfltS[ict][iz][iep][ik] << endl;
		}	}	}	}
	}


	// BBC
	// ------------------------------------------------------------------------
	TProfile *BBX[NordBBC][NsubBBC];
	TProfile *BBY[NordBBC][NsubBBC];
	for( int ith=0; ith<NordBBC; ith++ ){
		for( int iep=0; iep<NsubBBC; iep++ ){
			BBX[ith][iep] = (TProfile*)in->Get( Form("hBBCrecX_%d_%d",ith,iep) );
			BBY[ith][iep] = (TProfile*)in->Get( Form("hBBCrecY_%d_%d",ith,iep) );
			if( BBX[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get BBX!!"<< endl; return; }
			if( BBY[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get BBY!!"<< endl; return; }
		}
	}

	//cout << "BBC recentering parameter" << endl;
	for( int ith=0; ith<NordBBC; ith++ ){
		for( int ict=0; ict<NqvCent; ict++ ){
			for( int iz=0; iz<NqvZvtx; iz++ ){
				int id = 1+iz+ict*NqvZvtx;
				for( int iep=0; iep<NsubBBC; iep++ ){
					BBCXYmn[ith][ict][iz][iep][0] = BBX[ith][iep]->GetBinContent(id);
					BBCXYsg[ith][ict][iz][iep][0] = BBX[ith][iep]->GetBinError  (id);
					BBCXYmn[ith][ict][iz][iep][1] = BBY[ith][iep]->GetBinContent(id);
					BBCXYsg[ith][ict][iz][iep][1] = BBY[ith][iep]->GetBinError  (id);
					//cout <<ith<<" "<<ict<<" "<<iz<<" "<<iep<<" "<< BBCXYmn[ith][ict][iz][iep][0] <<" "<< BBCXYsg[ith][ict][iz][iep][0] <<" "<<  BBCXYmn[ith][ict][iz][iep][1] <<" "<< BBCXYsg[ith][ict][iz][iep][1] << endl;
	}	}	}	}
	//cout << endl;


	// Flattening parameters
	if( m_CalibMode==3 ){
		TProfile *fltBBCC[NordBBC][NsubBBC];
		TProfile *fltBBCS[NordBBC][NsubBBC];
		for( int ith=0; ith<NordBBC; ith++ ){
			for( int iep=0; iep<NsubBBC; iep++ ){
   				fltBBCC[ith][iep] = (TProfile*)in->Get( Form("hBBCfltC_%d_%d",ith,iep) );
   				fltBBCS[ith][iep] = (TProfile*)in->Get( Form("hBBCfltS_%d_%d",ith,iep) );
				if( fltBBCC[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltBBCC!!"<< endl; return; }
				if( fltBBCS[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltBBCS!!"<< endl; return; }
		}	}
		for( int ith=0; ith<NordBBC; ith++ ){
			for( int ict=0; ict<NqvCent; ict++ ){
				for( int iz=0; iz<NqvZvtx; iz++ ){
					for( int iep=0; iep<NsubBBC; iep++ ){
						for( int ik=0; ik<Nflt; ik++ ){
							int id = 1 + ik + iz*Nflt + ict*Nflt*NqvZvtx;
							BBCfltC[ith][ict][iz][iep][ik] = fltBBCC[ith][iep]->GetBinContent(id);
							BBCfltS[ith][ict][iz][iep][ik] = fltBBCS[ith][iep]->GetBinContent(id);
							//cout <<ith<<" "<<ict<<" "<<iz<<" "<<iep<<" "<< BBCfltC[ith][ict][iz][iep][ik] <<" "<< BBCfltS[ith][ict][iz][iep][ik] << endl;
		}	}	}	}	}
	}


	// TPC
	// ------------------------------------------------------------------------
	TProfile *hTPCX[NordTPC][NsubTPC];
	TProfile *hTPCY[NordTPC][NsubTPC];
	for( int ith=0; ith<NordTPC; ith++ ){
		for( int iep=0; iep<NsubTPC; iep++ ){
			hTPCX[ith][iep] = (TProfile*)in->Get( Form("hTPCrecX_%d_%d",ith,iep) );
			hTPCY[ith][iep] = (TProfile*)in->Get( Form("hTPCrecY_%d_%d",ith,iep) );
			if( hTPCX[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get hTPCX!!"<< endl; return; }
			if( hTPCY[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get hTPCY!!"<< endl; return; }
		}
	}

	//cout << "TPC recentering parameter" << endl;
	for( int ith=0; ith<NordTPC; ith++ ){
		for( int ict=0; ict<NqvCent; ict++ ){
			for( int iz=0; iz<NqvZvtx; iz++ ){
				int id = 1+iz+ict*NqvZvtx;
				for( int iep=0; iep<NsubTPC; iep++ ){
					TPCXYmn[ith][ict][iz][iep][0] = hTPCX[ith][iep]->GetBinContent(id);
					TPCXYsg[ith][ict][iz][iep][0] = hTPCX[ith][iep]->GetBinError  (id);
					TPCXYmn[ith][ict][iz][iep][1] = hTPCY[ith][iep]->GetBinContent(id);
					TPCXYsg[ith][ict][iz][iep][1] = hTPCY[ith][iep]->GetBinError  (id);
					//cout <<ith<<" "<<ict<<" "<<iz<<" "<<iep<<" "<< TPCXYmn[ith][ict][iz][iep][0] <<" "<< TPCXYsg[ith][ict][iz][iep][0] <<" "<<  TPCXYmn[ith][ict][iz][iep][1] <<" "<< TPCXYsg[ith][ict][iz][iep][1] << endl;
	}	}	}	}
	//cout << endl;

	if( m_CalibMode==3 ){
		TProfile *fltTPCC[NordTPC][NsubTPC];
		TProfile *fltTPCS[NordTPC][NsubTPC];
		for( int ith=0; ith<NordTPC; ith++ ){
			for( int iep=0; iep<NsubTPC; iep++ ){
				fltTPCC[ith][iep] = (TProfile*)in->Get( Form("hTPCfltC_%d_%d",ith,iep) );
				fltTPCS[ith][iep] = (TProfile*)in->Get( Form("hTPCfltS_%d_%d",ith,iep) );
				if( fltTPCC[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltTPCC!!"<< endl; return; }
				if( fltTPCS[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltTPCS!!"<< endl; return; }
		}	}
		for( int ith=0; ith<NordTPC; ith++ ){
			for( int ict=0; ict<NqvCent; ict++ ){
				for( int iz=0; iz<NqvZvtx; iz++ ){
					for( int iep=0; iep<NsubTPC; iep++ ){
						for( int ik=0; ik<Nflt; ik++ ){
							int id = 1 + ik + iz*Nflt +ict*Nflt*NqvZvtx;
							TPCfltC[ith][ict][iz][iep][ik] = fltTPCC[ith][iep]->GetBinContent(id);
							TPCfltS[ith][ict][iz][iep][ik] = fltTPCS[ith][iep]->GetBinContent(id);
							//cout <<ith<<" "<<ict<<" "<<iz<<" "<<iep<<" "<< TPCfltC[ith][ict][iz][iep][ik] <<" "<< TPCfltS[ith][ict][iz][iep][ik] << endl;
		}	}	}	}	}
	}

/*
	// EEMC
	// ------------------------------------------------------------------------
	TProfile *hEMCX[NordEMC];
	TProfile *hEMCY[NordEMC];
	for( int ith=0; ith<NordEMC; ith++ ){
		hEMCX[ith] = (TProfile*)in->Get( Form("hEMCrecX_%d",ith) );
		hEMCY[ith] = (TProfile*)in->Get( Form("hEMCrecY_%d",ith) );
		if( hEMCX[ith]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get hEMCX!!"<< endl; return; }
		if( hEMCY[ith]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get hEMCY!!"<< endl; return; }
	}

	for( int ith=0; ith<NordEMC; ith++ ){
		for( int ict=0; ict<NqvCent; ict++ ){
			for( int iz=0; iz<NqvZvtx; iz++ ){
				int id = 1+iz+ict*NqvZvtx;
				EMCXYmn[ith][ict][iz][0] = hEMCX[ith]->GetBinContent(id);
				EMCXYsg[ith][ict][iz][0] = hEMCX[ith]->GetBinError  (id);
				EMCXYmn[ith][ict][iz][1] = hEMCY[ith]->GetBinContent(id);
				EMCXYsg[ith][ict][iz][1] = hEMCY[ith]->GetBinError  (id);
				//if( ict==0 && iz==5 ) cout <<ith<<" "<<ict<<" "<<iz<<" "<< EMCXYmn[ith][ict][iz][0] <<" "<< EMCXYsg[ith][ict][iz][0] <<" "<<  EMCXYmn[ith][ict][iz][1] <<" "<< EMCXYsg[ith][ict][iz][1] << endl;
	}	}	}
	//cout << endl;

	if( m_CalibMode==3 ){
		TProfile *fltEMCC[NordEMC];
		TProfile *fltEMCS[NordEMC];
		for( int ith=0; ith<NordEMC; ith++ ){
			fltEMCC[ith] = (TProfile*)in->Get( Form("hEMCfltC_%d",ith) );
			fltEMCS[ith] = (TProfile*)in->Get( Form("hEMCfltS_%d",ith) );
			if( fltEMCC[ith]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltEMCC!!"<< endl; return; }
			if( fltEMCS[ith]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltEMCS!!"<< endl; return; }
		}

		for( int ith=0; ith<NordEMC; ith++ ){
			for( int ict=0; ict<NqvCent; ict++ ){
				for( int iz=0; iz<NqvZvtx; iz++ ){
					for( int ik=0; ik<Nflt; ik++ ){
						int id = 1 + ik + iz*Nflt +ict*Nflt*NqvZvtx;
						EMCfltC[ith][ict][iz][ik] = fltEMCC[ith]->GetBinContent(id);
						EMCfltS[ith][ict][iz][ik] = fltEMCS[ith]->GetBinContent(id);
						//if( ict==0 && iz==0 && ik<2 ) cout <<ith<<" "<<ict<<" "<<iz<<" "<<ik<<" "<< EMCfltC[ith][ict][iz][ik] <<" "<< EMCfltS[ith][ict][iz][ik] << endl;
		}	}	}	}
	}
*/
	// EPD
	// ------------------------------------------------------------------------
	TProfile *hEPDX[NordEPD][NsubEPD];
	TProfile *hEPDY[NordEPD][NsubEPD];
	for( int ith=0; ith<NordEPD; ith++ ){
		for( int iep=0; iep<NsubEPD; iep++ ){
			hEPDX[ith][iep] = (TProfile*)in->Get( Form("hEPDrecX_%d_%d",ith,iep) );
			hEPDY[ith][iep] = (TProfile*)in->Get( Form("hEPDrecY_%d_%d",ith,iep) );
			if( hEPDX[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get hEPDX!!"<< endl; return; }
			if( hEPDY[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get hEPDY!!"<< endl; return; }
		}
	}

	for( int ith=0; ith<NordEPD; ith++ ){
		for( int ict=0; ict<NqvCent; ict++ ){
			for( int iz=0; iz<NqvZvtx; iz++ ){
				int id = 1+iz+ict*NqvZvtx;
				for( int iep=0; iep<NsubEPD; iep++ ){
					EPDXYmn[ith][ict][iz][iep][0] = hEPDX[ith][iep]->GetBinContent(id);
					EPDXYsg[ith][ict][iz][iep][0] = hEPDX[ith][iep]->GetBinError  (id);
					EPDXYmn[ith][ict][iz][iep][1] = hEPDY[ith][iep]->GetBinContent(id);
					EPDXYsg[ith][ict][iz][iep][1] = hEPDY[ith][iep]->GetBinError  (id);
					//if( ict==0 && iz==5 ) cout <<ith<<" "<<ict<<" "<<iz<<" "<< EPDXYmn[ith][ict][iz][0] <<" "<< EPDXYsg[ith][ict][iz][0] <<" "<<  EPDXYmn[ith][ict][iz][1] <<" "<< EPDXYsg[ith][ict][iz][1] << endl;
	}	}	}	}
	//cout << endl;

	if( m_CalibMode==3 ){
		TProfile *fltEPDC[NordEPD][NsubEPD];
		TProfile *fltEPDS[NordEPD][NsubEPD];
		for( int ith=0; ith<NordEPD; ith++ ){
			for( int iep=0; iep<NsubEPD; iep++ ){
				fltEPDC[ith][iep] = (TProfile*)in->Get( Form("hEPDfltC_%d_%d",ith,iep) );
				fltEPDS[ith][iep] = (TProfile*)in->Get( Form("hEPDfltS_%d_%d",ith,iep) );
				if( fltEPDC[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltEPDC!!"<< endl; return; }
				if( fltEPDS[ith][iep]==NULL ){ cout <<"QVMaker::LoadEPCalbParam(): Could not get fltEPDS!!"<< endl; return; }
		}	}

		for( int ith=0; ith<NordEPD; ith++ ){
			for( int ict=0; ict<NqvCent; ict++ ){
				for( int iz=0; iz<NqvZvtx; iz++ ){
					for( int iep=0; iep<NsubTPC; iep++ ){
						for( int ik=0; ik<Nflt; ik++ ){
							int id = 1 + ik + iz*Nflt +ict*Nflt*NqvZvtx;
							EPDfltC[ith][ict][iz][iep][ik] = fltEPDC[ith][iep]->GetBinContent(id);
							EPDfltS[ith][ict][iz][iep][ik] = fltEPDS[ith][iep]->GetBinContent(id);
							//if( ict==0 && iz==0 && ik<2 ) cout <<ith<<" "<<ict<<" "<<iz<<" "<<ik<<" "<< EPDfltC[ith][ict][iz][ik] <<" "<< EPDfltS[ith][ict][iz][ik] << endl;
		}	}	}	}	}
	}


}

void QVMaker::setProfiles( TProfile *hrecX, TProfile *hrecY, TProfile *hfltC, TProfile *hfltS, TH1F *hEP[3] ){

	this->hrecX = hrecX;
	this->hrecY = hrecY;
	this->hfltC = hfltC;
	this->hfltS = hfltS;

	this->hEP[0] = hEP[0];
	this->hEP[1] = hEP[1];
	this->hEP[2] = hEP[2];
}
