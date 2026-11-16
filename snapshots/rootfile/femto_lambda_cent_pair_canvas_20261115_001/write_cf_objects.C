void write_cf_objects(const char* binsPath, const char* metadataPath, const char* rootPath);
#ifndef __CINT__
#include "../femto_lambda_ep_rapidity_cf_20261115_001/make_cf_plots.C"

void write_cf_objects(const char* binsPath, const char* metadataPath, const char* rootPath) {
  using namespace ep_rapidity_cf_plots;
  try {
    RootErrorGuard errors;
    gROOT->SetBatch(kTRUE); TH1::AddDirectory(kFALSE);
    Require(gSystem->AccessPathName(rootPath), "Refusing to overwrite CF objects");
    Selections selections; ReadSelections(binsPath,metadataPath,selections);
    RootOutput output(rootPath);
    for (Selections::const_iterator it=selections.begin();it!=selections.end();++it)
      WriteSelection(output.file,it->second);
    output.Close(); CheckRoot();
    std::cout << "WRITE_CF_OBJECTS_OK selections=" << selections.size() << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "WRITE_CF_OBJECTS_FAILED: " << error.what() << std::endl; gSystem->Exit(2);
  }
}
#endif
