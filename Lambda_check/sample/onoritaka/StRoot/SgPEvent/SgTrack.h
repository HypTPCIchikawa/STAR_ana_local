#ifndef SgTrack_hh
#define SgTrack_hh

//#include "TObject.h"
#include "StThreeVectorF.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h" // define kilogauss
#include "StMuDSTMaker/COMMON/StMuTrack.h"
//#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

//class SgTrack : public TObject {
class SgTrack {
	public:
		SgTrack();
		SgTrack(StMuTrack *mu);
		SgTrack(StPicoTrack *pico, StPicoEvent *peve, StPicoDst *pdst);
		~SgTrack();
		SgTrack& operator=(const SgTrack&);
		StPhysicalHelixD helix() const { return mHelix; }  // Returns inner helix (first measured point)
		//StThreeVectorF dcaGlobal(Int_t vtx_id=-1) const { return mGDca; }
		StThreeVector<double> momentum() const { return mMomVec; }
		StThreeVectorF tofPos() const { return mTofPos; }
		float dca   () const { return mGDca; }
		float beta  () const { return mTofBeta; }
		float tof   () const { return mTof; }
		float mom   () const { return mMom; }
		short charge() const { return mCH;  }
		short flg   () const { return mFlg; }

		void SetHelix(StPhysicalHelixD& h);
		//void SetGDca(StThreeVectorF& v);
		void SetTofBeta(float val);
		void set_flg(short val){ mFlg = val; }

	private:
		StPhysicalHelixD mHelix;
		StThreeVector<double> mMomVec;
		StThreeVectorF mTofPos;
		float mGDca;
		float mTofBeta;
		float mTof;
		float mMom;
		short mCH;
		short mFlg; // 0=likely BG,  1=likely lambda

	//ClassDef(SgTrack,0);
};
typedef std::vector<SgTrack> Trks;
typedef std::vector<Trks> TrksBuff;

#endif /* SgTrack_hh */
