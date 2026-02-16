#ifndef SgPEvent_hh
#define SgPEvent_hh

//#include "SgTrack/SgTrack.h"
#include "SgTrack.h"
#include "StThreeVectorF.hh"
#include <vector>

class SgPEvent {
	public:
		SgPEvent();
		~SgPEvent();
		SgPEvent& operator=(const SgPEvent&);
		//SgPEvent(const SgPEvent&);

		void set_pVtx(StThreeVectorF& vtx){ mpVtx = vtx; }
		void set_BField(double val){ mBField = val; }
		void set_Cent(int val){ mCent = val; };
		void set_Pid (int val){ mPid = val; };
		void set_Ch  (int val){ mCh  = val; }
		void set_Psi(float val[3]){ 
			for(int k=0; k<3; k++) mPsiSP[k] = val[k];
		}
		void set_Psi2(float val[3]){ 
			for(int k=0; k<3; k++) mPsi2T[k] = val[k];
	   	}
		void set_Psi2B(float val[3]){ 
			for(int k=0; k<3; k++) mPsi2B[k] = val[k];
	   	}
		void set_Q1B(float val){ mQ1B = val; }
		void set_Q1E(float val){ mQ1E = val; }
		void set_Trgeff(float val){ mTrgEff = val; }
		void set_ResZDC(int i, float val){ mResZDC[i] = val; }
		void set_trks(Trks& ftrk){ trk = ftrk; }

		StThreeVectorF get_pVtx() const { return mpVtx; }
		double get_BField()       const { return mBField; }
		int get_Cent()            const { return mCent; }
		int get_Pid()             const { return mPid; }
		int get_Ch()              const { return mCh; }
		float get_Psi(int i)      const { return mPsiSP[i]; }
		float get_Psi2(int i)     const { return mPsi2T[i]; }
		float get_Psi2B(int i)    const { return mPsi2B[i]; }
		float get_Q1B()           const { return mQ1B; }
		float get_Q1E()           const { return mQ1E; }
		float get_Trgeff()        const { return mTrgEff; }
		float get_ResZDC(int i)   const { return mResZDC[i]; }
		Trks& get_trks(){ return trk; }

		void reserveTrkCap(int n);
		void AddTrk(SgTrack& sgl_trk);
		void Clear();

	private:
		StThreeVectorF mpVtx;
		double mBField;
		int    mCent;
		int    mPid; // 0=pion 1=kaon 2=proton
		int    mCh;
		float  mPsiSP[3]; // spectator planes
		float  mPsi2T[3];  // TPC 2nd-order event plane
		float  mPsi2B[3];  // BBC 2nd-order event plane
		float  mQ1B;    // BBC E+W Q1
		float  mQ1E;    // EEMC Q1
		float  mTrgEff;
		float  mResZDC[2];

		Trks trk;
};
typedef std::vector<SgPEvent> PEvents;
#endif
