#include "QV.h"

QV::QV() :  m_Det(dZDC), m_Sub(0), m_Ord(0), m_Cent(0), m_Zvtx(0) {

	ResetQ();
}

QV::~QV(){
}

void QV::set_Qv(float& fQx, float& fQy, float& fQw, DETID fDet, int& fSub, int& fOrd, int& fCent, int& fVz){

	m_Qx = fQx;
	m_Qy = fQy;
	m_Qw = fQw;
	m_Qabs = sqrt( m_Qx*m_Qx + m_Qy*m_Qy );

	m_Det  = fDet;
	m_Sub  = fSub;
	m_Ord  = fOrd;
	m_Cent = fCent;
	m_Zvtx = fVz;
}

void QV::set_Qv(float& fQx, float& fQy, float& fQw, float& fQmult, DETID fDet, int& fSub, int& fOrd, int& fCent, int& fVz){

	m_Qx    = fQx;
	m_Qy    = fQy;
	m_Qw    = fQw;
	m_Qmult = fQmult;
	m_Qabs = sqrt( m_Qx*m_Qx + m_Qy*m_Qy );

	m_Det  = fDet;
	m_Sub  = fSub;
	m_Ord  = fOrd;
	m_Cent = fCent;
	m_Zvtx = fVz;
}

void QV::ResetQ(){

	m_Qx     = -9999.0;
	m_Qy     = -9999.0;
	m_Qabs   = -9999.0;
	m_Qw     = -9999.0;
	m_Psi    = -9999.0;
	m_Qmult  = -9999.0;
}
