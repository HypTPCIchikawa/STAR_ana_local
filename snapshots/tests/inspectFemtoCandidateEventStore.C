#include "../include/FemtoCandidateEventStore.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

void inspectFemtoCandidateEventStore(const char* input,const char* report) {
 using namespace femto_candidate_store;
 try {
  if(!input||!*input||!report||!*report||!gSystem->AccessPathName(report))throw std::runtime_error("new explicit store/report paths required");
  TFile file(input,"READ");if(file.IsZombie())throw std::runtime_error("cannot open store");
  TObjString* layout=dynamic_cast<TObjString*>(file.Get(kLayoutName));
  TObjString* selection=dynamic_cast<TObjString*>(file.Get(kSelectionName));TObjString* mixing=dynamic_cast<TObjString*>(file.Get(kMixingName));
  TTree* tree=dynamic_cast<TTree*>(file.Get(kTreeName));if(!layout||!selection||!mixing||!tree)throw std::runtime_error("missing candidate store");
  std::vector<std::string> species;std::istringstream lines(layout->GetString().Data());std::string line;
  while(std::getline(lines,line)){std::istringstream words(line);std::string key,value;words>>key>>value;if(key=="species")species.push_back(value);}
  const std::string selected=selection->GetString().Data(),mixed=mixing->GetString().Data(),uuid=file.GetUUID().AsString();
  const Long64_t declared=tree->GetEntries(),compressed=tree->GetZipBytes(),uncompressed=tree->GetTotBytes();file.Close();
  Reader reader;if(!reader.Init(input,species,selected,mixed,true,true))throw std::runtime_error(reader.LastError());
  Event e;FemtoCandidateStore physical;FemtoNuclearPidPairQaStore qa;
  Long64_t rows=0,empty=0,qaOnly=0,pool=0,otherSpeciesOnly=0,originalEp=0,storedEp=0,qaRows=0,qaSelected=0;std::map<std::string,Long64_t> counts;
  std::set<int> runs;int previous=-1;double zmin=1.e99,zmax=-1.e99,bmin=1.e99,bmax=-1.e99;
  while(reader.Next(e,physical,qa)) {
   if(reader.NewFile())previous=-1;
   if(e.inputIndex<=previous)throw std::runtime_error("input order is not strictly increasing within original stream");previous=e.inputIndex;++rows;runs.insert(e.runId);
   zmin=std::min(zmin,e.pvZ);zmax=std::max(zmax,e.pvZ);bmin=std::min(bmin,e.bField);bmax=std::max(bmax,e.bField);
   if(e.eventPlaneValid)++originalEp;
   if(e.eventPlane2Valid){if(e.eventPlane2Tracks<=0||e.eventPlane2<0||e.eventPlane2>=TMath::Pi()||e.q2x*e.q2x+e.q2y*e.q2y<=0)throw std::runtime_error("invalid auxiliary TPC Q2");++storedEp;}
   Long64_t n=0,paired=0,nqa=0;
   for(FemtoCandidateStore::const_iterator i=physical.begin();i!=physical.end();++i){counts[i->first]+=i->second.size();n+=i->second.size();if(i->first=="lambda"||i->first=="he4"||i->first=="he3")paired+=i->second.size();}
   for(FemtoNuclearPidPairQaStore::const_iterator i=qa.begin();i!=qa.end();++i)for(size_t j=0;j<i->second.size();++j){++qaRows;++nqa;if(i->second[j].selected)++qaSelected;}
   if(n)++pool;else{++empty;if(nqa)++qaOnly;}if(n&&paired==0)++otherSpeciesOnly;
  }
  if(!reader.Finished()||rows!=declared)throw std::runtime_error(reader.LastError().empty()?"cache row count mismatch":reader.LastError());
  std::ofstream out(report);if(!out)throw std::runtime_error("cannot create cache inspection report");out<<std::setprecision(17);
  out<<"{\"phase\":\"passed\",\"rows\":"<<rows<<",\"originalOrderVerified\":true,\"physicalCandidates\":{";bool comma=false;
  for(std::map<std::string,Long64_t>::const_iterator i=counts.begin();i!=counts.end();++i){if(comma)out<<',';comma=true;out<<'\"'<<i->first<<"\":"<<i->second;}
  out<<"},\"eventsWithPhysicalCandidates\":"<<pool<<",\"eventsWithoutPhysicalCandidates\":"<<empty<<",\"qaOnlyEvents\":"<<qaOnly<<",\"eventsWithOnlyUnpairedDT\":"<<otherSpeciesOnly<<",\"qaCandidates\":"<<qaRows<<",\"qaSelectedCandidates\":"<<qaSelected<<",\"distinctRuns\":"<<runs.size()<<",\"originalValidEventPlanes\":"<<originalEp<<",\"auxiliaryValidEventPlanes\":"<<storedEp<<",\"auxiliaryInvalidEventPlanes\":"<<rows-storedEp<<",\"vertexZRangeCm\":["<<zmin<<','<<zmax<<"],\"fieldRangeKGauss\":["<<bmin<<','<<bmax<<"],\"treeZipBytes\":"<<compressed<<",\"treeUncompressedBytes\":"<<uncompressed<<",\"sourceFileUUID\":\""<<uuid<<"\"}\n";
  if(!out)throw std::runtime_error("cache inspection write failed");std::cout<<"PASS inspected "<<rows<<" cached accepted events, valid auxiliary Q2 "<<storedEp<<std::endl;
 } catch(const std::exception& e){std::cerr<<"CANDIDATE_CACHE_INSPECTION_FAILED "<<e.what()<<std::endl;gSystem->Exit(2);}
}
