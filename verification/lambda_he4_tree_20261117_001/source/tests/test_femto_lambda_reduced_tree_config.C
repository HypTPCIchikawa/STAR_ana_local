int test_femto_lambda_reduced_tree_config(const char*,const char*,const char*);
#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "include/FemtoLambdaReducedTree.h"
#include "include/FemtoLambdaReducedTreeConfig.h"
#include "TString.h"
#include <iostream>
#include <stdexcept>

int test_femto_lambda_reduced_tree_config(const char* compactMain,const char* fullMain,const char* fullPilot) {
  try {
    if(!compactMain||!*compactMain||!fullMain||!*fullMain||!fullPilot||!*fullPilot)
      throw std::runtime_error("explicit compact/full mainconfs and completed full-control pilot are required");
    femto_lambda_tree::Settings compact,full;compact.Load(compactMain);full.Load(fullMain);
    if(!compact.enabled||!full.enabled||compact.storePrePidControl||!full.storePrePidControl)
      throw std::runtime_error("compact/full profiles do not specify their intended control-bank capability");
    if(compact.selectionProvenance!=full.selectionProvenance+femto_lambda_tree::kPrePidControlOmittedMarker)
      throw std::runtime_error("compact fingerprint must differ by the explicit omitted-control marker only");
    const std::string fixture=std::string(gSystem->WorkingDirectory())+"/tmp/reduced-config-missing-"+TString::Format("%d",gSystem->GetPid()).Data();
    if(gSystem->mkdir(fixture.c_str(),kTRUE)!=0)throw std::runtime_error("cannot create missing-key fixture directory");
    const YAML::Node originalMain=YAML::LoadFile(compactMain);
    const YAML::Node originalFeature=YAML::LoadFile(femto_candidate_store::ConfigReference(compactMain,originalMain["reducedTree"].as<std::string>()));
    YAML::Node badFeature,badMain;
    for(YAML::const_iterator i=originalFeature.begin();i!=originalFeature.end();++i)
      if(i->first.as<std::string>()!="storePrePidControl")badFeature[i->first.as<std::string>()]=i->second;
    badMain["reducedTree"]=fixture+"/feature.yaml";
    {std::ofstream output((fixture+"/feature.yaml").c_str());output<<YAML::Dump(badFeature);if(!output)throw std::runtime_error("missing-key feature write failed");}
    {std::ofstream output((fixture+"/main.yaml").c_str());output<<YAML::Dump(badMain);if(!output)throw std::runtime_error("missing-key main write failed");}
    bool rejectedMissing=false;
    try {femto_lambda_tree::Settings missing;missing.Load(fixture+"/main.yaml");}
    catch(const std::exception& error){rejectedMissing=std::string(error.what()).find("missing storePrePidControl")!=std::string::npos;}
    if(!rejectedMissing)throw std::runtime_error("missing storePrePidControl silently accepted");
    femto_lambda_tree::Reader historical;
    if(!historical.Init(fullPilot,full.species,full.selectionProvenance,full.mixingProvenance,true))
      throw std::runtime_error("existing full-control pilot compatibility failed: "+historical.LastError());
    femto_lambda_tree::Reader mismatched;
    if(mismatched.Init(fullPilot,compact.species,compact.selectionProvenance,compact.mixingProvenance,true))
      throw std::runtime_error("full-control payload was accepted as a compact-envelope input");
    std::cout<<"PASS reduced-tree capability profiles: mandatory flag, distinct compact envelope, unchanged full-control pilot fingerprint ("
             <<historical.EntriesTotal()<<" entries)"<<std::endl;
    return 0;
  } catch(const std::exception& error) {std::cerr<<"FAIL reduced-tree capability config: "<<error.what()<<std::endl;return 1;}
}
#endif
