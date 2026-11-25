#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include <iostream>

void CheckSuspectRoot(const char* first,const char* second) {
  const char* inputs[]={first,second};
  for(int i=0;i<2;++i) {
    TFile f(inputs[i],"READ");
    std::cout<<"EARLY_PICO_KEY_AUDIT path="<<inputs[i]<<" zombie="<<f.IsZombie()
             <<" recovered="<<f.TestBit(TFile::kRecovered)<<" fileBytes="<<f.GetSize()
             <<" logicalEnd="<<f.GetEND()<<std::endl;
    if(f.IsZombie()) continue;
    TIter keys(f.GetListOfKeys());
    while(TObject* key=keys()) {
      TKey* k=dynamic_cast<TKey*>(key);
      if(k) std::cout<<"key="<<k->GetName()<<" class="<<k->GetClassName()<<std::endl;
    }
    std::cout<<"picoKey="<<(f.GetKey("PicoDst")!=0)
             <<" picoTree="<<(dynamic_cast<TTree*>(f.Get("PicoDst"))!=0)<<std::endl;
  }
}
