#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
void RunProbePhiStarHelix(const char* settings,const char* output) {
 const char* base="/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_audit_20261111_001/source/ProbePhiStarHelix.C";
 if(gSystem->Load("libStarClassLibrary.so")<0) {std::cerr<<"Cannot load StarClassLibrary"<<std::endl;gSystem->Exit(2);return;}
 gSystem->AddIncludePath(" -I/star/nfs4/AFS/star/packages/SL24y/StRoot/StarClassLibrary");
 gSystem->AddLinkedLibs(" -L/star/nfs4/AFS/star/packages/SL24y/.sl73_x8664_gcc485/lib -lStarClassLibrary");
 if(gROOT->ProcessLine(Form(".L %s+",base))!=0) {gSystem->Exit(2);return;}
 gROOT->ProcessLine(Form("ProbePhiStarHelix(\"%s\",\"%s\")",settings,output));
}
