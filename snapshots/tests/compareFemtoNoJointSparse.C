// Compare a normal writer with an otherwise identical writer whose joint PID
// Sparse is explicitly disabled. Both inputs are opened read-only.
#include "compareFemtoCandidateReplay.C"
#include "../include/FemtoCandidateEventStore.h"
#include "TClass.h"
#include "TROOT.h"
#include <cstring>

namespace FemtoNoJointComparison {
using FemtoReplayComparison::Require;
using FemtoReplayComparison::Number;
using FemtoReplayComparison::Quote;

void HistogramKeys(TDirectory& directory, const std::string& prefix,
                   std::set<std::string>& dense, std::set<std::string>& sparse) {
  TIter next(directory.GetListOfKeys());
  TKey* key = 0;
  while ((key = static_cast<TKey*>(next()))) {
    const std::string name = key->GetName();
    const std::string path = prefix + name;
    const std::string className = key->GetClassName();
    if (className == "TDirectoryFile" || className == "TDirectory") {
      TDirectory* child = directory.GetDirectory(name.c_str());
      Require(child != 0, "cannot open directory: " + path);
      HistogramKeys(*child, path + "/", dense, sparse);
      continue;
    }
    TClass* type = gROOT->GetClass(className.c_str());
    if (type && type->InheritsFrom(TH1::Class()))
      Require(dense.insert(path).second, "duplicate dense histogram cycle: " + path);
    if ((type && type->InheritsFrom(THnBase::Class())) ||
        className.find("THnSparse") != std::string::npos)
      Require(sparse.insert(path).second, "duplicate multidimensional histogram cycle: " + path);
  }
}

void AxisPresentation(const TAxis& a, const TAxis& b, const std::string& path) {
  FemtoReplayComparison::Axis(a, b, path);
  Require(std::string(a.GetTitle()) == b.GetTitle(), "axis title: " + path);
  Require(a.GetFirst() == b.GetFirst() && a.GetLast() == b.GetLast(),
          "axis display range: " + path);
  for (Int_t i = 1; i <= a.GetNbins(); ++i)
    Require(std::string(a.GetBinLabel(i)) == b.GetBinLabel(i),
            "axis bin label: " + path + " bin " + Number(i));
}

template<class T> bool SameBits(const T& a, const T& b) {
  return std::memcmp(&a, &b, sizeof(T)) == 0;
}

template<class T> void VectorBits(const std::vector<T>& a, const std::vector<T>& b,
                                  const std::string& name, Long64_t row) {
  Require(a.size() == b.size(), "cache vector length: " + name + " row " + Number(row));
  if (!a.empty())
    Require(std::memcmp(&a[0], &b[0], a.size() * sizeof(T)) == 0,
            "cache vector payload: " + name + " row " + Number(row));
}

void EventBits(const femto_candidate_store::Event& a,
               const femto_candidate_store::Event& b, Long64_t row) {
#define FCNS_EVENT_FIELD(field) \
  Require(SameBits(a.field, b.field), "cache event field " #field " row " + Number(row))
  FCNS_EVENT_FIELD(inputIndex);
  FCNS_EVENT_FIELD(runId);
  FCNS_EVENT_FIELD(eventId);
  FCNS_EVENT_FIELD(cent9);
  FCNS_EVENT_FIELD(cent16);
  FCNS_EVENT_FIELD(refMult);
  FCNS_EVENT_FIELD(fxtMult);
  FCNS_EVENT_FIELD(nTracks);
  FCNS_EVENT_FIELD(nBTOFMatch);
  FCNS_EVENT_FIELD(eventPlane2Tracks);
  FCNS_EVENT_FIELD(pvX);
  FCNS_EVENT_FIELD(pvY);
  FCNS_EVENT_FIELD(pvZ);
  FCNS_EVENT_FIELD(bField);
  FCNS_EVENT_FIELD(refMultCorr);
  FCNS_EVENT_FIELD(centWeight);
  FCNS_EVENT_FIELD(psi2);
  FCNS_EVENT_FIELD(eventPlane2);
  FCNS_EVENT_FIELD(q2x);
  FCNS_EVENT_FIELD(q2y);
  FCNS_EVENT_FIELD(eventPlaneValid);
  FCNS_EVENT_FIELD(eventPlane2Valid);
  FCNS_EVENT_FIELD(epdPsi1);
  FCNS_EVENT_FIELD(epdRawQx);
  FCNS_EVENT_FIELD(epdRawQy);
  FCNS_EVENT_FIELD(epdRecenteredQx);
  FCNS_EVENT_FIELD(epdRecenteredQy);
  FCNS_EVENT_FIELD(epdWeight);
  FCNS_EVENT_FIELD(epdStatus);
  FCNS_EVENT_FIELD(epdBin);
  FCNS_EVENT_FIELD(epdRunIndex);
  FCNS_EVENT_FIELD(epdCalibrationCent);
  FCNS_EVENT_FIELD(epdHits);
  FCNS_EVENT_FIELD(epdValid);
#undef FCNS_EVENT_FIELD
}

std::set<std::string> BranchNames(TTree& tree) {
  std::set<std::string> names;
  TIter next(tree.GetListOfBranches());
  TBranch* branch = 0;
  while ((branch = static_cast<TBranch*>(next())))
    Require(names.insert(branch->GetName()).second, "duplicate cache branch");
  return names;
}

void SameStoreMetadata(TFile& a, TFile& b) {
  const char* names[] = {
    femto_candidate_store::kLayoutName,
    femto_candidate_store::kSelectionName,
    femto_candidate_store::kMixingName
  };
  for (unsigned i = 0; i < sizeof(names)/sizeof(names[0]); ++i) {
    TObjString* left = dynamic_cast<TObjString*>(a.Get(names[i]));
    TObjString* right = dynamic_cast<TObjString*>(b.Get(names[i]));
    Require(left && right && left->GetString() == right->GetString(),
            std::string("cache metadata mismatch: ") + names[i]);
  }
  // UUIDs must be different for independently produced files and each must
  // describe its own original stream, rather than a copied or hadd-ed tree.
  TObjString* left = dynamic_cast<TObjString*>(a.Get(femto_candidate_store::kFileUuidName));
  TObjString* right = dynamic_cast<TObjString*>(b.Get(femto_candidate_store::kFileUuidName));
  Require(left && right, "both original-stream UUID markers are required");
  Require(left->GetString() == a.GetUUID().AsString() &&
          right->GetString() == b.GetUUID().AsString(), "cache file UUID marker mismatch");
  Require(left->GetString() != right->GetString(), "independent producers must have distinct UUIDs");
}

Long64_t Cache(TTree& a, TTree& b, Long64_t expectedRows,
               Long64_t& candidates, Long64_t& floatValues) {
  using namespace femto_candidate_store;
  Require(a.GetEntries() == expectedRows && b.GetEntries() == expectedRows,
          "cache row count does not match explicit expected count");
  Require(BranchNames(a) == BranchNames(b), "cache branch sets must match exactly");
  Event eventA, eventB;
  Buffers bufferA, bufferB;
  std::string error;
  const bool epd = a.GetBranch("epdPsi1") != 0;
  Require(Bind(&a, eventA, bufferA, false, error, epd), "reference cache binding: " + error);
  error.clear();
  Require(Bind(&b, eventB, bufferB, false, error, epd), "new cache binding: " + error);
  for (Long64_t row = 0; row < expectedRows; ++row) {
    Require(a.GetEntry(row) > 0 && b.GetEntry(row) > 0,
            "cache read failed at row " + Number(row));
    EventBits(eventA, eventB, row);
    VectorBits(*bufferA.species, *bufferB.species, "species", row);
    VectorBits(*bufferA.collection, *bufferB.collection, "collection", row);
    VectorBits(*bufferA.qaSelected, *bufferB.qaSelected, "qaSelected", row);
    VectorBits(*bufferA.physicalPresent, *bufferB.physicalPresent, "physicalPresent", row);
    VectorBits(*bufferA.qaPresent, *bufferB.qaPresent, "qaPresent", row);
    VectorBits(*bufferA.floats, *bufferB.floats, "Float payload", row);
    VectorBits(*bufferA.ints, *bufferB.ints, "Int payload", row);
    VectorBits(*bufferA.quality, *bufferB.quality, "quality flags", row);
    candidates += static_cast<Long64_t>(bufferA.species->size());
    floatValues += static_cast<Long64_t>(bufferA.floats->size());
  }
  a.ResetBranchAddresses();
  b.ResetBranchAddresses();
  return expectedRows;
}

void Completed(TFile& file) {
  TNamed* status = dynamic_cast<TNamed*>(file.Get("FemtoLambdaRunStatus"));
  Require(status && std::string(status->GetTitle()) == "completed",
          "both producers must finish with status completed");
}

Double_t SignalEntries(TFile& file, const char* name, Double_t expected) {
  TH1* histogram = dynamic_cast<TH1*>(file.Get(name));
  Require(histogram && histogram->GetEntries() == expected && expected > 0,
          std::string("nonzero expected He4 signal entries required: ") + name);
  return histogram->GetEntries();
}
} // namespace FemtoNoJointComparison

void compareFemtoNoJointSparse(const char* reference, const char* output,
                              const char* report, Long64_t expectedCacheRows,
                              Double_t expectedHe4SE, Double_t expectedHe4ME) {
  using namespace FemtoNoJointComparison;
  try {
    Require(reference && *reference && output && *output && report && *report,
            "explicit input and new report paths required");
    Require(expectedCacheRows > 0 && expectedHe4SE > 0 && expectedHe4ME > 0,
            "explicit positive cache and exercised pair counts required");
    Require(gSystem->AccessPathName(report), "refusing existing comparison report");
    TFile a(reference, "READ"), b(output, "READ");
    Require(!a.IsZombie() && !b.IsZombie() && !a.TestBit(TFile::kRecovered) &&
            !b.TestBit(TFile::kRecovered), "valid closed ROOT inputs required");
    Completed(a);
    Completed(b);
    std::set<std::string> denseA, denseB, sparseA, sparseB;
    HistogramKeys(a, "", denseA, sparseA);
    HistogramKeys(b, "", denseB, sparseB);
    Require(!denseA.empty() && denseA == denseB, "all dense histogram key sets must match");
    Require(sparseA.size() == 1 && sparseA.count("hNuclearPidPairJoint") == 1,
            "reference must contain exactly the joint PID Sparse");
    Require(sparseB.empty() && !b.GetKey("hNuclearPidPairJoint"),
            "new output must contain no multidimensional Sparse histogram");
    Require(!b.GetKey("NuclearPidPairJointMetadata"),
            "disabled Sparse must not retain stale joint metadata");
    SignalEntries(a, "hKstarSE_lambda_he4_signal", expectedHe4SE);
    SignalEntries(a, "hKstarME_lambda_he4_signal", expectedHe4ME);
    const Double_t se = SignalEntries(b, "hKstarSE_lambda_he4_signal", expectedHe4SE);
    const Double_t me = SignalEntries(b, "hKstarME_lambda_he4_signal", expectedHe4ME);
    Long64_t denseCells = 0;
    for (std::set<std::string>::const_iterator i = denseA.begin(); i != denseA.end(); ++i) {
      TH1* left = dynamic_cast<TH1*>(a.Get(i->c_str()));
      TH1* right = dynamic_cast<TH1*>(b.Get(i->c_str()));
      Require(left && right && std::string(left->ClassName()) == right->ClassName(),
              "dense class mismatch: " + *i);
      Require(std::string(left->GetTitle()) == right->GetTitle(), "histogram title: " + *i);
      AxisPresentation(*left->GetXaxis(), *right->GetXaxis(), *i + " X");
      if (left->GetDimension() > 1)
        AxisPresentation(*left->GetYaxis(), *right->GetYaxis(), *i + " Y");
      if (left->GetDimension() > 2)
        AxisPresentation(*left->GetZaxis(), *right->GetZaxis(), *i + " Z");
      denseCells += FemtoReplayComparison::Dense(*left, *right, *i);
      delete left;
      delete right;
    }
    TTree* treeA = dynamic_cast<TTree*>(a.Get(femto_candidate_store::kTreeName));
    TTree* treeB = dynamic_cast<TTree*>(b.Get(femto_candidate_store::kTreeName));
    Require(treeA && treeB, "both complete candidate trees are required");
    SameStoreMetadata(a, b);
    Long64_t candidates = 0, floatValues = 0;
    const Long64_t rows = Cache(*treeA, *treeB, expectedCacheRows, candidates, floatValues);
    std::ofstream result(report);
    Require(bool(result), "cannot create comparison report");
    result << std::setprecision(17)
      << "{\n  \"phase\": \"passed\",\n"
      << "  \"reference\": " << Quote(reference) << ",\n"
      << "  \"output\": " << Quote(output) << ",\n"
      << "  \"denseHistograms\": " << denseA.size() << ",\n"
      << "  \"denseCellsIncludingFlows\": " << denseCells << ",\n"
      << "  \"exactClassesAxesTitlesLabelsContentsErrorsSumw2EntriesStatistics\": true,\n"
      << "  \"sparseHistogramsReference\": " << sparseA.size() << ",\n"
      << "  \"sparseHistogramsOutput\": 0,\n"
      << "  \"jointMetadataAbsent\": true,\n"
      << "  \"bothProducerStatusesCompleted\": true,\n"
      << "  \"he4SignalSEEntries\": " << se << ",\n"
      << "  \"he4SignalMEEntries\": " << me << ",\n"
      << "  \"cacheRows\": " << rows << ",\n"
      << "  \"cacheCandidatesIncludingQa\": " << candidates << ",\n"
      << "  \"cacheFloatValues\": " << floatValues << ",\n"
      << "  \"cacheEventAndCandidatePayloadBitwiseEqual\": true,\n"
      << "  \"cacheSelectionMixingAndLayoutEqual\": true,\n"
      << "  \"originalStreamUuidsVerified\": true\n}\n";
    Require(bool(result), "comparison report write failed");
    std::cout << "PASS no-joint Sparse: " << denseA.size() << " dense histograms, "
              << denseCells << " cells; " << rows << " bitwise-equal cache rows; He4 SE="
              << se << " ME=" << me << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "NO_JOINT_SPARSE_COMPARISON_FAILED " << error.what() << std::endl;
    gSystem->Exit(2);
  }
}
