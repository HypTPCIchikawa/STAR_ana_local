#include <TFile.h>
#include <TH2.h>
#include <TKey.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
void audit_he4_y(const char* listPath,const char* outputPath) {
  std::ifstream paths(listPath); if(!paths)throw std::runtime_error("Missing explicit input list");
  std::vector<double> edges,counts;double entries=0;std::vector<double> sourceEntries;std::vector<std::string> names;
  std::string path;
  while(std::getline(paths,path)){
    if(path.empty()||path[0]=='#')continue;
    TFile input(path.c_str(),"READ");if(input.IsZombie())throw std::runtime_error("Unreadable input");
    TH2* h=dynamic_cast<TH2*>(input.Get("hNucleus_PtVsYLab_he4"));if(!h)throw std::runtime_error("Missing He4 single-particle QA");
    if(edges.empty()){
      for(int b=1;b<=h->GetNbinsX()+1;++b)edges.push_back(h->GetXaxis()->GetBinLowEdge(b));
      counts.assign(h->GetNbinsX()+2,0.);
    }
    if(counts.size()!=static_cast<size_t>(h->GetNbinsX()+2))throw std::runtime_error("Axis mismatch");
    for(int b=1;b<=h->GetNbinsX()+1;++b)if(std::fabs(edges[b-1]-h->GetXaxis()->GetBinLowEdge(b))>1e-12)throw std::runtime_error("Axis edges changed");
    for(int b=0;b<=h->GetNbinsX()+1;++b)for(int y=0;y<=h->GetNbinsY()+1;++y)counts[b]+=h->GetBinContent(b,y);
    entries+=h->GetEntries();sourceEntries.push_back(h->GetEntries());names.push_back(path);
  }
  if(names.empty())throw std::runtime_error("Empty input list");
  std::ofstream out(outputPath);if(!out)throw std::runtime_error("Cannot create output");
  out<<std::setprecision(17)<<"{\n\"histogram\":\"hNucleus_PtVsYLab_he4\",\"xDefinition\":\"selected physical He4 P4.Rapidity() laboratory frame\",\"population\":\"fixed producer single-particle QA; all accepted centralities and EP validity states; not a pair-weighted/cent9-filtered distribution\",\"sourceCount\":"<<names.size()<<",\"entries\":"<<entries<<",\"edges\":[";
  for(size_t i=0;i<edges.size();++i){if(i)out<<',';out<<edges[i];}out<<"],\"countsIncludingPtFlow\":[";
  for(size_t i=0;i<counts.size();++i){if(i)out<<',';out<<counts[i];}out<<"],\"sources\":[";
  for(size_t i=0;i<names.size();++i){if(i)out<<',';out<<"{\"path\":\""<<names[i]<<"\",\"entries\":"<<sourceEntries[i]<<'}';}out<<"],\"treeEntriesRead\":0}\n";
  out.close();if(!out)throw std::runtime_error("Write failed");
  std::cout<<"HE4_Y_AUDIT_PASS sources="<<names.size()<<" entries="<<entries<<std::endl;
}
