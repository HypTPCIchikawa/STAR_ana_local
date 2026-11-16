#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TNamed.h"
#include "TSystem.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// Compare the two already-completed cached replays. No reconstruction or replay
// runs here. The input fixture holds two real events in the same original pool:
// reference event IDs are 1,2; two-file IDs are 1,1 before the Maker remaps them.
namespace FemtoIdentityReplayComparison {

void Require(bool ok, const std::string& message)
{
  if (!ok) throw std::runtime_error(message);
}

TH1* Get(TFile& file, const char* name)
{
  TH1* histogram = dynamic_cast<TH1*>(file.Get(name));
  Require(histogram != 0, std::string("missing histogram ") + name);
  return histogram;
}

Long64_t CellsIncludingFlows(const TH1& histogram)
{
  Long64_t count = histogram.GetNbinsX() + 2;
  if (histogram.GetDimension() >= 2) count *= histogram.GetNbinsY() + 2;
  if (histogram.GetDimension() >= 3) count *= histogram.GetNbinsZ() + 2;
  return count;
}

void CompareAxis(const TAxis& a, const TAxis& b)
{
  Require(a.GetNbins() == b.GetNbins(), "histogram axis bin count differs");
  for (Int_t bin = 1; bin <= a.GetNbins() + 1; ++bin) {
    Require(a.GetBinLowEdge(bin) == b.GetBinLowEdge(bin),
            "histogram axis bin edge differs");
  }
}

Long64_t Compare(TH1& a, TH1& b)
{
  Require(std::string(a.ClassName()) == b.ClassName() &&
          a.GetDimension() == b.GetDimension(),
          std::string("histogram type/dimension differs: ") + a.GetName());
  CompareAxis(*a.GetXaxis(), *b.GetXaxis());
  if (a.GetDimension() >= 2) CompareAxis(*a.GetYaxis(), *b.GetYaxis());
  if (a.GetDimension() >= 3) CompareAxis(*a.GetZaxis(), *b.GetZaxis());
  const Long64_t cells = CellsIncludingFlows(a);
  Require(cells == CellsIncludingFlows(b), "histogram cell count differs");
  Require(a.GetEntries() == b.GetEntries(),
          std::string("histogram Entries differ: ") + a.GetName());
  for (Long64_t bin = 0; bin < cells; ++bin) {
    Require(a.GetBinContent(bin) == b.GetBinContent(bin) &&
            a.GetBinError(bin) == b.GetBinError(bin),
            std::string("histogram content/error differs including flows: ") +
              a.GetName());
  }
  Double_t statsA[13] = {0};
  Double_t statsB[13] = {0};
  a.GetStats(statsA);
  b.GetStats(statsB);
  for (Int_t statistic = 0; statistic < 13; ++statistic) {
    Require(statsA[statistic] == statsB[statistic],
            std::string("histogram statistics differ: ") + a.GetName());
  }
  return cells;
}

void Run(const char* reference, const char* split, const char* output,
         Long64_t expected, Int_t statusOverlap, Int_t statusEligible,
         Int_t statusAttempted, Int_t statusFilled)
{
  Require(reference && *reference && split && *split && output && *output &&
          expected > 0, "explicit outputs and positive expected ME required");
  Require(gSystem->AccessPathName(output),
          "refuse existing identity comparison JSON");
  TFile a(reference, "READ");
  TFile b(split, "READ");
  Require(!a.IsZombie() && !b.IsZombie(),
          "cannot open identity replay outputs");
  TNamed* statusA = dynamic_cast<TNamed*>(a.Get("FemtoLambdaRunStatus"));
  TNamed* statusB = dynamic_cast<TNamed*>(b.Get("FemtoLambdaRunStatus"));
  Require(statusA && statusB &&
          std::string(statusA->GetTitle()) == "completed" &&
          std::string(statusB->GetTitle()) == "completed",
          "identity replay producer not completed");

  const char* names[] = {
    "hKstarSE_lambda_he4_signal",
    "hKstarSEVsCent_lambda_he4_signal",
    "hKstarME_lambda_he4_signal",
    "hKstarMEVsCent_lambda_he4_signal",
    "hMixSamplerQA"
  };
  Long64_t cells = 0;
  for (unsigned histogram = 0;
       histogram < sizeof(names) / sizeof(names[0]); ++histogram) {
    cells += Compare(*Get(a, names[histogram]), *Get(b, names[histogram]));
  }

  TH1* meA = Get(a, "hKstarME_lambda_he4_signal");
  TH1* meB = Get(b, "hKstarME_lambda_he4_signal");
  Require(meA->GetEntries() == expected && meB->GetEntries() == expected,
          "positive collision-witness ME signal count not reproduced");
  TH2* qaA = dynamic_cast<TH2*>(Get(a, "hMixSamplerQA"));
  TH2* qaB = dynamic_cast<TH2*>(Get(b, "hMixSamplerQA"));
  Require(qaA && qaB, "sampler QA is not two-dimensional");

  // Existing sampler axis convention is x=status, y=channel. Use Double_t
  // coordinates to avoid ROOT5's FindBin(const char*) overload ambiguity at 0.
  const Int_t channel = qaA->GetYaxis()->FindBin(Double_t(0));
  const Int_t overlap = qaA->GetXaxis()->FindBin(Double_t(statusOverlap));
  const Int_t eligible = qaA->GetXaxis()->FindBin(Double_t(statusEligible));
  const Int_t attempted = qaA->GetXaxis()->FindBin(Double_t(statusAttempted));
  const Int_t filled = qaA->GetXaxis()->FindBin(Double_t(statusFilled));
  Require(channel >= 1 && channel <= qaA->GetNbinsY() &&
          overlap >= 1 && overlap <= qaA->GetNbinsX() &&
          eligible >= 1 && eligible <= qaA->GetNbinsX() &&
          attempted >= 1 && attempted <= qaA->GetNbinsX() &&
          filled >= 1 && filled <= qaA->GetNbinsX(),
          "sampler status/channel coordinates fall outside the booked axes");
  Require(qaA->GetBinContent(overlap, channel) == 0 &&
          qaB->GetBinContent(overlap, channel) == 0,
          "false cross-event overlap veto remains");
  Require(qaA->GetBinContent(eligible, channel) > 0 &&
          qaA->GetBinContent(eligible, channel) ==
            qaA->GetBinContent(attempted, channel) &&
          qaA->GetBinContent(eligible, channel) ==
            qaA->GetBinContent(filled, channel),
          "reference sampler did not fill every eligible pair");
  Require(qaB->GetBinContent(eligible, channel) ==
            qaB->GetBinContent(attempted, channel) &&
          qaB->GetBinContent(eligible, channel) ==
            qaB->GetBinContent(filled, channel),
          "split sampler did not fill every eligible pair");

  std::ofstream report(output);
  Require(report.good(), "cannot create identity comparison JSON");
  report << std::setprecision(17)
    << "{\n  \"phase\":\"passed\",\n"
    << "  \"denseHistogramsCompared\":5,\n"
    << "  \"cellsIncludingFlows\":" << cells << ",\n"
    << "  \"exactAxesContentErrorsEntriesStatistics\":true,\n"
    << "  \"producerStatusReference\":\"completed\",\n"
    << "  \"producerStatusSplit\":\"completed\",\n"
    << "  \"expectedHe4MeSignalEntries\":" << expected << ",\n"
    << "  \"referenceHe4MeSignalEntries\":" << meA->GetEntries() << ",\n"
    << "  \"splitHe4MeSignalEntries\":" << meB->GetEntries() << ",\n"
    << "  \"referenceHe4SignalSEEntries\":"
    << Get(a, "hKstarSE_lambda_he4_signal")->GetEntries() << ",\n"
    << "  \"splitHe4SignalSEEntries\":"
    << Get(b, "hKstarSE_lambda_he4_signal")->GetEntries() << ",\n"
    << "  \"samplerXAxis\":\"status\",\n"
    << "  \"samplerYAxis\":\"channel\",\n"
    << "  \"samplerXAxisBins\":" << qaA->GetNbinsX() << ",\n"
    << "  \"samplerXAxisMin\":" << qaA->GetXaxis()->GetXmin() << ",\n"
    << "  \"samplerXAxisMax\":" << qaA->GetXaxis()->GetXmax() << ",\n"
    << "  \"samplerYAxisBins\":" << qaA->GetNbinsY() << ",\n"
    << "  \"samplerYAxisMin\":" << qaA->GetYaxis()->GetXmin() << ",\n"
    << "  \"samplerYAxisMax\":" << qaA->GetYaxis()->GetXmax() << ",\n"
    << "  \"samplerHe4ChannelCoordinate\":0,\n"
    << "  \"samplerHe4ChannelBin\":" << channel << ",\n"
    << "  \"samplerOverlapStatusCoordinate\":" << statusOverlap << ",\n"
    << "  \"samplerOverlapStatusBin\":" << overlap << ",\n"
    << "  \"samplerOverlapVetoReference\":"
    << qaA->GetBinContent(overlap, channel) << ",\n"
    << "  \"samplerOverlapVetoSplit\":"
    << qaB->GetBinContent(overlap, channel) << ",\n"
    << "  \"samplerEligibleStatusCoordinate\":" << statusEligible << ",\n"
    << "  \"samplerAttemptedStatusCoordinate\":" << statusAttempted << ",\n"
    << "  \"samplerFilledStatusCoordinate\":" << statusFilled << ",\n"
    << "  \"samplerEligibleReference\":"
    << qaA->GetBinContent(eligible, channel) << ",\n"
    << "  \"samplerAttemptedReference\":"
    << qaA->GetBinContent(attempted, channel) << ",\n"
    << "  \"samplerFilledReference\":"
    << qaA->GetBinContent(filled, channel) << ",\n"
    << "  \"samplerEligibleSplit\":"
    << qaB->GetBinContent(eligible, channel) << ",\n"
    << "  \"samplerAttemptedSplit\":"
    << qaB->GetBinContent(attempted, channel) << ",\n"
    << "  \"samplerFilledSplit\":"
    << qaB->GetBinContent(filled, channel) << ",\n"
    << "  \"samplerEligibleAttemptedFilledAllMatch\":true,\n"
    << "  \"all20FloatFieldsFixedInInputFixture\":true,\n"
    << "  \"noPicoDstOrKfReconstruction\":true\n}\n";
  report.close();
  Require(!report.fail(), "identity comparison JSON close failed");
  std::cout << "FEMTO_CANDIDATE_IDENTITY_REPLAY_PASS ME="
            << meA->GetEntries() << " cells=" << cells << std::endl;
}

} // namespace FemtoIdentityReplayComparison

void compareFemtoCandidateIdentityReplay(
    const char* reference, const char* split, const char* output,
    Long64_t expectedMeSignal, Int_t statusOverlap, Int_t statusEligible,
    Int_t statusAttempted, Int_t statusFilled)
{
  try {
    FemtoIdentityReplayComparison::Run(reference, split, output,
      expectedMeSignal, statusOverlap, statusEligible, statusAttempted,
      statusFilled);
  } catch (const std::exception& error) {
    std::cerr << "FEMTO_CANDIDATE_IDENTITY_REPLAY_FAILED: "
              << error.what() << std::endl;
    gSystem->Exit(2);
  }
}
