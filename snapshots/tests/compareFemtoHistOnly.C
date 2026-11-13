// Compare a candidate-store writer with a histogram-only producer using the
// same physical selection. The inputs are read-only; only a new report is made.
#include "compareFemtoNoJointSparse.C"
#include "TParameter.h"

namespace FemtoHistOnlyComparison {
using FemtoReplayComparison::Require;
using FemtoReplayComparison::Number;
using FemtoReplayComparison::Quote;

Long64_t Counter(TFile& file, const char* name, Long64_t expected) {
  TParameter<Long64_t>* value =
      dynamic_cast<TParameter<Long64_t>*>(file.Get(name));
  Require(value && value->GetVal() == expected,
          std::string("explicit expected event counter required: ") + name);
  return value->GetVal();
}

void NoCandidateStore(TFile& file) {
  const char* keys[] = {
    femto_candidate_store::kTreeName,
    femto_candidate_store::kLayoutName,
    femto_candidate_store::kSelectionName,
    femto_candidate_store::kMixingName,
    femto_candidate_store::kFileUuidName,
    "FemtoCandidateEventStoreMode",
    "FemtoCandidateEventStoreSettings",
    "FemtoCandidateEventStoreInput",
    "FemtoCandidateReplayScope"
  };
  for (unsigned i = 0; i < sizeof(keys) / sizeof(keys[0]); ++i)
    Require(!file.GetKey(keys[i]),
            std::string("histogram-only output retains candidate-store key: ") + keys[i]);
}
} // namespace FemtoHistOnlyComparison

void compareFemtoHistOnly(const char* reference, const char* output,
                         const char* report, Long64_t expectedDense,
                         Long64_t expectedInput, Long64_t expectedAccepted,
                         Double_t expectedHe4SE, Double_t expectedHe4ME) {
  using namespace FemtoHistOnlyComparison;
  try {
    Require(reference && *reference && output && *output && report && *report,
            "explicit reference, output and new report paths required");
    Require(expectedDense > 0 && expectedInput > 0 && expectedAccepted > 0 &&
            expectedAccepted <= expectedInput && expectedHe4SE > 0 && expectedHe4ME > 0,
            "explicit positive histogram, event and exercised pair counts required");
    Require(gSystem->AccessPathName(report), "refusing existing comparison report");
    TFile a(reference, "READ"), b(output, "READ");
    Require(!a.IsZombie() && !b.IsZombie() && !a.TestBit(TFile::kRecovered) &&
            !b.TestBit(TFile::kRecovered), "valid closed ROOT inputs required");
    FemtoNoJointComparison::Completed(a);
    FemtoNoJointComparison::Completed(b);
    Counter(a, "inputEvents", expectedInput);
    Counter(b, "inputEvents", expectedInput);
    Counter(a, "reconstructedEvents", expectedAccepted);
    Counter(b, "reconstructedEvents", expectedAccepted);

    TTree* cache = dynamic_cast<TTree*>(a.Get(femto_candidate_store::kTreeName));
    Require(cache && cache->GetEntries() == expectedAccepted,
            "reference must contain all expected accepted candidate-store rows");
    TNamed* mode = dynamic_cast<TNamed*>(a.Get("FemtoCandidateEventStoreMode"));
    Require(mode && std::string(mode->GetTitle()) == "write",
            "reference must be a completed candidate-store writer");
    NoCandidateStore(b);

    std::set<std::string> denseA, denseB, sparseA, sparseB;
    FemtoNoJointComparison::HistogramKeys(a, "", denseA, sparseA);
    FemtoNoJointComparison::HistogramKeys(b, "", denseB, sparseB);
    Require(static_cast<Long64_t>(denseA.size()) == expectedDense && denseA == denseB,
            "all recursive dense histogram keys must match the explicit expected count");
    Require(sparseA.empty() && sparseB.empty(),
            "both configurations must disable every multidimensional Sparse histogram");
    Require(!a.GetKey("hNuclearPidPairJoint") && !b.GetKey("hNuclearPidPairJoint") &&
            !a.GetKey("NuclearPidPairJointMetadata") &&
            !b.GetKey("NuclearPidPairJointMetadata"),
            "joint Sparse or stale joint metadata must not be present");
    FemtoNoJointComparison::SignalEntries(a, "hKstarSE_lambda_he4_signal", expectedHe4SE);
    FemtoNoJointComparison::SignalEntries(a, "hKstarME_lambda_he4_signal", expectedHe4ME);
    const Double_t se = FemtoNoJointComparison::SignalEntries(
        b, "hKstarSE_lambda_he4_signal", expectedHe4SE);
    const Double_t me = FemtoNoJointComparison::SignalEntries(
        b, "hKstarME_lambda_he4_signal", expectedHe4ME);

    Long64_t cells = 0;
    for (std::set<std::string>::const_iterator i = denseA.begin(); i != denseA.end(); ++i) {
      TH1* left = dynamic_cast<TH1*>(a.Get(i->c_str()));
      TH1* right = dynamic_cast<TH1*>(b.Get(i->c_str()));
      Require(left && right && std::string(left->ClassName()) == right->ClassName(),
              "dense class mismatch: " + *i);
      Require(std::string(left->GetTitle()) == right->GetTitle(),
              "histogram title mismatch: " + *i);
      FemtoNoJointComparison::AxisPresentation(
          *left->GetXaxis(), *right->GetXaxis(), *i + " X");
      if (left->GetDimension() > 1)
        FemtoNoJointComparison::AxisPresentation(
            *left->GetYaxis(), *right->GetYaxis(), *i + " Y");
      if (left->GetDimension() > 2)
        FemtoNoJointComparison::AxisPresentation(
            *left->GetZaxis(), *right->GetZaxis(), *i + " Z");
      cells += FemtoReplayComparison::Dense(*left, *right, *i);
      delete left;
      delete right;
    }

    std::ofstream result(report);
    Require(bool(result), "cannot create histogram-only comparison report");
    result << std::setprecision(17)
      << "{\n  \"phase\": \"passed\",\n"
      << "  \"reference\": " << Quote(reference) << ",\n"
      << "  \"output\": " << Quote(output) << ",\n"
      << "  \"denseHistograms\": " << denseA.size() << ",\n"
      << "  \"denseCellsIncludingFlows\": " << cells << ",\n"
      << "  \"exactClassesAxesTitlesLabelsContentsErrorsSumw2EntriesStatistics\": true,\n"
      << "  \"bothProducerStatusesCompleted\": true,\n"
      << "  \"inputEvents\": " << expectedInput << ",\n"
      << "  \"reconstructedEvents\": " << expectedAccepted << ",\n"
      << "  \"he4SignalSEEntries\": " << se << ",\n"
      << "  \"he4SignalMEEntries\": " << me << ",\n"
      << "  \"sparseHistogramsReference\": 0,\n"
      << "  \"sparseHistogramsOutput\": 0,\n"
      << "  \"jointMetadataAbsent\": true,\n"
      << "  \"referenceCandidateStoreRows\": " << cache->GetEntries() << ",\n"
      << "  \"candidateStoreRowsOutput\": 0,\n"
      << "  \"candidateStoreTreeAbsent\": true,\n"
      << "  \"candidateStoreModeAbsent\": true,\n"
      << "  \"referenceRootBytes\": " << a.GetEND() << ",\n"
      << "  \"histOnlyRootBytes\": " << b.GetEND() << ",\n"
      << "  \"outputScopeChange\": \"Candidate event store disabled; dense histograms unchanged\"\n}\n";
    Require(bool(result), "comparison report write failed");
    std::cout << "PASS histogram-only closure: " << denseA.size() << " dense histograms, "
              << cells << " exact cells; candidate store absent; He4 SE=" << se
              << " ME=" << me << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "HIST_ONLY_COMPARISON_FAILED " << error.what() << std::endl;
    gSystem->Exit(2);
  }
}
