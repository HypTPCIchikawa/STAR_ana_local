#include "SgPEvent.h"

SgPEvent::SgPEvent(){
	Clear();
}

SgPEvent::~SgPEvent(){
}

SgPEvent& SgPEvent::operator=(const SgPEvent& ev){

	if( this != &ev ){
		this->mpVtx   = ev.get_pVtx();
		this->mBField = ev.get_BField();
		this->mCent   = ev.get_Cent();
		this->mPid    = ev.get_Pid();
		this->mCh     = ev.get_Ch();
		for(int k=0; k<3; k++) this->mPsiSP[k] = ev.get_Psi(k);
		for(int k=0; k<3; k++) this->mPsi2B[k] = ev.get_Psi2B(k);
		for(int k=0; k<3; k++) this->mPsi2T[k] = ev.get_Psi2(k);
		this->mQ1B    = ev.get_Q1B();
		this->mQ1E    = ev.get_Q1E();
		this->mTrgEff = ev.get_Trgeff();
		for(int k=0; k<2; k++) this->mResZDC[k]  = ev.get_ResZDC(k);

		this->trk     = ev.trk;
		//this->trk     = ev.get_trks();
	}

	return *this;
}

void SgPEvent::reserveTrkCap(int n){
	trk.reserve(n);
}

void SgPEvent::AddTrk(SgTrack& sgl_trk){
	trk.push_back( sgl_trk );
}

void SgPEvent::Clear(){

	mpVtx   = StThreeVectorF(-9999,-9999,-9999);
	mBField = -9999.;
	mCent   = -1;
	mPid    = -1;
	mCh     = -9999;
	for( int iep=0; iep<3; iep++ ){
		mPsiSP[iep] = -9999.;
		mPsi2B[iep] = -9999.;
		mPsi2T[iep] = -9999.;
	}
	mQ1B    = -9999.;
	mQ1E    = -9999.;
	mTrgEff = 1.0;
	for( int k=0; k<2; k++ ) mResZDC[k] = 1.0;
	trk.clear();
}
