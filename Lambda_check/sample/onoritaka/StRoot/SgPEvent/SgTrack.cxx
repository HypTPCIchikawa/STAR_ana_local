#include "SgTrack.h"

// tof
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "phys_constants.h"

//ClassImp(SgTrack);

SgTrack::SgTrack(){
}

//SgTrack::SgTrack(StMuTrack *mu){
//	mHelix   = mu->helix();
//	mGDca    = mu->dcaGlobal().mag();
//	mTofPos  = mu->btofPidTraits().position();
//	mTofBeta = mu->btofPidTraits().beta();
//	mTof     = mu->btofPidTraits().timeOfFlight();
//	mMom     = mu->p().mag();
//	mMomVec  = mu->p();
//	mCH      = mu->charge();
//	mFlg     = 0;
//}

SgTrack::SgTrack(StPicoTrack *pico, StPicoEvent *peve, StPicoDst *pdst){


	// run11
	//StPhysicalHelixD helix( pico->gMom(), pico->origin(), BField*kilogauss, pico->charge() );
	//mGDca    = pico->dca();
	//mHelix   = helix;
	//mMom     = pico->gMom().mag();
	//mMomVec  = pico->gMom();
	//mTofPos  = pico->btofHisPos();
	//mTofBeta = pico->btofBeta();
	//mTof     = pico->btof();
	  
	// run14
	//mHelix   = pico->helix();
	//mGDca    = fabs( mHelix.geometricSignedDistance(peve->primaryVertex()) );
	//mMom     = pico->gMom( peve->primaryVertex(), peve->bField() ).mag();
	//mMomVec  = pico->gMom( peve->primaryVertex(), peve->bField() );
	//mCH      = pico->charge();
	//mFlg     = 0;

	// run16
	//mHelix   = pico->helix( peve->bField() );
	//mGDca    = pico->dca().mag();
	//mMom     = pico->gMom().mag();
	//mMomVec  = pico->gMom();
	//mCH      = pico->charge();
	//mFlg     = 0;

	// run18
	//mGDca      = pico->dca().mag();
	mGDca      = pico->gDCA( peve->primaryVertex() ).Mag();

	mMom       = pico->gMom().Mag();
	mCH        = pico->charge();
	mFlg       = 0;

	//mHelix     = pico->helix( peve->bField() );
	StThreeVectorF gmom( pico->gMom().X(), pico->gMom().Y(), pico->gMom().Z() );
	StThreeVectorF org( pico->origin().X(), pico->origin().Y(), pico->origin().Z() );
	StPhysicalHelixD helix( gmom, org, peve->bField()*kilogauss, static_cast<float>(pico->charge()) );
	mHelix = helix;

	//mMomVec    = pico->gMom();
	mMomVec    = gmom;

	// run14 test
	//StDcaGeometry dcaG = pico->dcaGeometry();
	//mMom    = dcaG.momentum().mag();
	//mMomVec = dcaG.momentum();

	//StThreeVectorF vtx = peve->primaryVertex();
	TVector3 vtx = peve->primaryVertex();

	int index = pico->bTofPidTraitsIndex();
	StPicoBTofPidTraits *btof;
   	if( index>=0 ) btof = pdst->btofPidTraits(index);

	if( index>=0 && btof ){
		mTof     = btof->btof();
		mTofBeta = btof->btofBeta();

		//mTofPos  = btof->btofHitPos();
		TVector3 tofPos3d  = btof->btofHitPos();
		StThreeVectorF tofPos( tofPos3d.X(), tofPos3d.Y(), tofPos3d.Z() );
		mTofPos = tofPos;

		if( mTofBeta<1e-4 ){

			StThreeVectorF vtxSt( vtx.X(), vtx.Y(), vtx.Z() );

			//float L   = tofPathLength(&vtx, &mTofPos, pico->helix(peve->bField()).curvature());
			float L   = tofPathLength(&vtxSt, &tofPos, pico->helix(peve->bField()).curvature());
			if( mTof>0 ) mTofBeta = L / (mTof * (C_C_LIGHT / 1.e9));
			else         mTofBeta = -1.0; //std::numeric_limits<float>::quiet_NaN();
		}
	}else{
		mTof     = -1.0;
		mTofPos  = StThreeVector<double>(-9999.,-9999.,-9999.);
		mTofBeta = -1.0;
	}
}

SgTrack::~SgTrack(){
}

SgTrack& SgTrack::operator=(const SgTrack& ftrk){

	if( this != &ftrk ){
		this->mHelix   = ftrk.mHelix;
		this->mGDca    = ftrk.mGDca;
		this->mTofPos  = ftrk.mTofPos;
		this->mTofBeta = ftrk.mTofBeta;
		this->mTof     = ftrk.mTof;
		this->mMom     = ftrk.mMom;
		this->mMomVec  = ftrk.mMomVec;
		this->mCH      = ftrk.mCH;
		this->mFlg     = ftrk.mFlg;
	}
	return *this;
}

void SgTrack::SetHelix(StPhysicalHelixD& h){
	mHelix = h;
}

//void SgTrack::SetGDca(StThreeVectorF& v){
//	mGDca = v;
//}

void SgTrack::SetTofBeta(float fbeta){
	mTofBeta = fbeta;
}
