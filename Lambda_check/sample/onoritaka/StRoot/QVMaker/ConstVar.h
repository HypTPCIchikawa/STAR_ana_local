#ifndef _CONSTVAR_H_
#define _CONSTVAR_H_

static const int NsubSMD = 3;
static const int NsubBBC = 6;
static const int NsubTPC = 13; //16;
static const int NsubEPD = 3;

static const int NordSMD = 2;
static const int NordBBC = 2;
static const int NordTPC = 5;
static const int NordEMC = 5;
static const int NordEPD = 2;

static const int NqvCent = 16; // 5% step upto 80%
static const int NqvZvtx = 10; // 6cm step within (-30,30) for Cu+Au or ...
//static const int NqvZvtx = 6; // 2cm step within (-6,6) for run14 AuAu200

enum DETID { dZDC, dBBC, dTPC, dEEMC, dEPD };
#endif
