#ifndef _QV_H_
#define _QV_H_

#include "ConstVar.h"
#include <cmath>

class QV {
	public:
		QV();
		~QV();

		void ResetQ();

		void set_Qx  (float val){ m_Qx   = val;  }
		void set_Qy  (float val){ m_Qy   = val;  }
		void set_Qabs(float val){ m_Qabs = val;  }
		void set_Qw  (float val){ m_Qw   = val;  }
		void set_Psi (float val){ m_Psi  = val; }
		void set_Ord (int val)  { m_Ord  = val; }
		void set_Qv  (float& fQx, float& fQy, float& fQw, DETID fDet, int& fSub, int& fOrd, int& fCent, int& fVz);
		void set_Qv (float& fQx, float& fQy, float& fQw, float& fQmult, DETID fDet, int& fSub, int& fOrd, int& fCent, int& fVz);

		float get_Qx    (){ return m_Qx;    }
		float get_Qy    (){ return m_Qy;    }
		float get_Qabs  (){ return m_Qabs;    }
		float get_Qw    (){ return m_Qw;    }
		float get_Qmult (){ return m_Qmult; }
		float get_Psi   (){ return m_Psi;   }
		DETID get_Det   (){ return m_Det;   }
		int   get_Sub   (){ return m_Sub;   }
		int   get_Ord   (){ return m_Ord;   }
		int   get_Cent  (){ return m_Cent;  }
		int   get_Zvtx  (){ return m_Zvtx;  }

	private:
		DETID m_Det; // 0:ZDC-SMD, 1:BBC, 2:TPC, 3:EEMC
		int m_Sub;
		int m_Ord;
		int m_Cent;
		int m_Zvtx;

		float m_Qx;
		float m_Qy;
		float m_Qabs;
		float m_Qw;
		float m_Psi;
		float m_Qmult;
};

#endif
