#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

void InventoryEarlyEpd(const char* cacheList,const char* auditPath,const char* expectedManifest) {
  try {
    std::ifstream files(cacheList);
    std::ofstream audit(auditPath);
    std::string path;
    int count=0;Long64_t total=0;
    while(std::getline(files,path)) {
      std::unique_ptr<TFile> source(TFile::Open(path.c_str(),"READ"));
      if(!source.get() || source->IsZombie() || source->TestBit(TFile::kRecovered))
        throw std::runtime_error("bad output "+path);
      TNamed* inventory=dynamic_cast<TNamed*>(source->Get("EpdCalibrationInputInventory"));
      TNamed* sha=dynamic_cast<TNamed*>(source->Get("EpdCalibrationManifestSHA256"));
      TNamed* id=dynamic_cast<TNamed*>(source->Get("EpdCalibrationJob"));
      TTree* raw=dynamic_cast<TTree*>(source->Get("epdRaw"));
      if(!inventory || !sha || !id || !raw || std::string(sha->GetTitle())!=expectedManifest ||
         !raw->GetBranch("q") || !raw->GetBranch("fold") || !raw->GetBranch("run"))
        throw std::runtime_error("bad collector provenance "+path);
      std::istringstream rows(inventory->GetTitle());
      std::string line;
      Long64_t kept=0;
      while(std::getline(rows,line)) {
        std::istringstream row(line);
        std::string leaf,extra; int run,fold; Long64_t nominal,entries,seen,valid;
        if(!(row>>leaf>>run>>fold>>nominal>>entries>>seen>>valid) || (row>>extra) ||
           entries<=0 || seen!=entries || valid<0 || valid>seen)
           throw std::runtime_error("bad per-file inventory "+path);
        kept+=valid;
        audit<<id->GetTitle()<<'\t'<<line<<'\n';
      }
      if(kept!=raw->GetEntries()) throw std::runtime_error("cache count mismatch "+path);
      total+=kept;++count;
    }
    if(!files.eof() || !count || !audit.good()) throw std::runtime_error("empty or failed inventory");
    std::cout<<"EPD_CACHE_AUDIT_PASS roots="<<count<<" rawEntries="<<total<<std::endl;
  } catch(const std::exception& ex) {
    std::cerr<<"EPD_CACHE_AUDIT_ERROR "<<ex.what()<<std::endl;
    gSystem->Exit(1);
  }
}
