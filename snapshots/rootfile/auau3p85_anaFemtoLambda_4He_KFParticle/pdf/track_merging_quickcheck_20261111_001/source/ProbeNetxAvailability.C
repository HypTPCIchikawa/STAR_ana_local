#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
void ProbeNetxAvailability(const char* url,const char* output){
 const char* base="/star/nfs4/AFS/star/ROOT/5.34.38/.sl73_x8664_gcc485/rootdeb/lib";
 gSystem->AddDynamicPath(base);
 const int loaded=gSystem->Load("/star/nfs4/AFS/star/ROOT/5.34.38/.sl73_x8664_gcc485/rootdeb/lib/libNetx.so");
 std::ofstream out(output);out<<"{\"rootVersion\":\""<<gROOT->GetVersion()<<"\",\"plugin\":\"same-version ROOT5.34.38 SL73 libNetx.so\",\"loadReturnCode\":"<<loaded;
 if(loaded<0){out<<",\"status\":\"load_failed\",\"picoOpened\":false}\n";out.close();std::cout<<"NETX_ISOLATED_PROBE_LOAD_FAILED"<<std::endl;return;}
 gPluginMgr->AddHandler("TFile","^[x]?root:","TXNetFile","/star/nfs4/AFS/star/ROOT/5.34.38/.sl73_x8664_gcc485/rootdeb/lib/libNetx.so","TXNetFile(const char*,Option_t*,const char*,Int_t,Int_t,Bool_t)");
 TFile* f=TFile::Open(url,"READ");TTree* t=(f&&!f->IsZombie())?dynamic_cast<TTree*>(f->Get("PicoDst")):0;
 if(!t){out<<",\"status\":\"remote_header_open_failed\",\"picoOpened\":false}\n";out.close();std::cout<<"NETX_ISOLATED_PROBE_OPEN_FAILED"<<std::endl;if(f)delete f;return;}
 out<<",\"status\":\"passed\",\"picoOpened\":true,\"treeEntries\":"<<t->GetEntries()<<"}\n";out.close();delete f;std::cout<<"NETX_ISOLATED_PROBE_PASS"<<std::endl;
}
