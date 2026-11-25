// Read-only QA for merged/new Lambda-nucleus ROOT output. No Phi mass model is reused.
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TObjString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"
#include "ConfigManager.h"
#include "cuts/FemtoConfig.h"
#include <iostream>
#include <string>
#include <vector>

namespace FemtoLambdaQa {
TH2* SumMass(TFile& file, const TString& pattern, const char* name, Int_t first, Int_t last) {
  TH2* sum = 0;
  for (Int_t cent = first; cent <= last; ++cent) {
    TH2* histogram = dynamic_cast<TH2*>(file.Get(TString::Format(pattern, cent)));
    if (!histogram) { delete sum; return 0; }
    if (!sum) { sum = static_cast<TH2*>(histogram->Clone(name)); sum->SetDirectory(0); }
    else sum->Add(histogram);
  }
  return sum;
}
void Message(const char* text) {
  TLatex label; label.SetNDC(); label.SetTextSize(0.035); label.DrawLatex(0.12, 0.5, text);
}

Int_t CheckAcceptanceAndSplit(TFile& file, const TString& channel,
                            const TString& species, TH2* acceptance[4]) {
  // Old ROOT output has none of the additions. A partial contract must not look successful.
  const TString names[] = {
    "hLambda_PtVsYLab_signal", "hLambdaProton_PtVsYLab_signal",
    "hLambdaPion_PtVsYLab_signal", "hNucleus_PtVsYLab_" + species,
    "hKstarSE_" + channel + "_signalLow", "hKstarSEVsCent_" + channel + "_signalLow",
    "hKstarME_" + channel + "_signalLow", "hKstarMEVsCent_" + channel + "_signalLow",
    "hKstarSE_" + channel + "_signalHigh", "hKstarSEVsCent_" + channel + "_signalHigh",
    "hKstarME_" + channel + "_signalHigh", "hKstarMEVsCent_" + channel + "_signalHigh"
  };
  const Int_t dimensions[] = {2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2};
  Int_t present = 0;
  for (Int_t i = 0; i < 12; ++i) {
    TObject* object = file.Get(names[i]);
    if (!object) continue;
    ++present;
    TH1* histogram = dynamic_cast<TH1*>(object);
    if (!histogram || histogram->GetDimension() != dimensions[i]) {
      std::cerr << "ERROR: wrong histogram type/dimension for " << names[i] << std::endl;
      return -1;
    }
    if (i < 4) acceptance[i] = dynamic_cast<TH2*>(histogram);
  }
  if (present == 0) {
    std::cerr << "WARNING: new acceptance/split-signal histograms absent; rerun analysis "
              << "to create them. Drawing the original four QA pages only." << std::endl;
    return 0;
  }
  if (present != 12) {
    std::cerr << "ERROR: incomplete acceptance/split-signal histogram contract ("
              << present << "/12); rerun analysis with matching Maker and YAML."
              << std::endl;
    return -1;
  }
  return 1;
}
void DrawAcceptance(TCanvas& canvas, TH2* acceptance[4], const char* outputPdf) {
  canvas.Clear(); canvas.Divide(2, 2);
  for (Int_t i = 0; i < 4; ++i) {
    canvas.cd(i + 1);
    acceptance[i]->SetStats(kFALSE);
    acceptance[i]->Draw("COLZ");
  }
  canvas.Print(outputPdf);
}

Bool_t SameAxis(const TAxis& first, const TAxis& second) {
  if (first.GetNbins() != second.GetNbins()) return kFALSE;
  for (Int_t bin = 1; bin <= first.GetNbins() + 1; ++bin)
    if (first.GetBinLowEdge(bin) != second.GetBinLowEdge(bin)) return kFALSE;
  return kTRUE;
}

Int_t CheckPairMt(TFile& file, const TString& shortName) {
  // Check all stored centralities, not only the range selected for the QA display.
  // Old output has no mT maps; a partial or incompatible schema is an error.
  const char* suffixes[] = {"", "_SBPos", "_SBNeg", "_signalLow", "_signalHigh"};
  Int_t present = 0;
  TH2* reference = 0;
  for (Int_t mixed = 0; mixed < 2; ++mixed) {
    const TString prefix = mixed ? "mix/hKstarMt_Mixed_" : "true/hKstarMt_";
    for (Int_t region = 0; region < 5; ++region) {
      for (Int_t cent = 0; cent < 9; ++cent) {
        const TString name = prefix + shortName + suffixes[region]
            + TString::Format("_CentBin%d", cent);
        TObject* object = file.Get(name);
        if (!object) continue;
        ++present;
        TH2* histogram = dynamic_cast<TH2*>(object);
        if (!histogram || histogram->GetDimension() != 2) {
          std::cerr << "ERROR: wrong histogram type/dimension for " << name << std::endl;
          return -1;
        }
        if (!reference) reference = histogram;
        else if (!SameAxis(*reference->GetXaxis(), *histogram->GetXaxis())
                 || !SameAxis(*reference->GetYaxis(), *histogram->GetYaxis())) {
          std::cerr << "ERROR: incompatible kstar-mT axes for " << name << std::endl;
          return -1;
        }
      }
    }
  }
  if (present == 0) {
    std::cerr << "WARNING: kstar-mT histograms absent; rerun analysis to create them. "
              << "Skipping the five kstar-mT QA pages." << std::endl;
    return 0;
  }
  if (present != 90) {
    std::cerr << "ERROR: incomplete kstar-mT histogram contract (" << present
              << "/90); rerun analysis with matching Maker and YAML." << std::endl;
    return -1;
  }
  return 1;
}

Int_t DrawPairMt(TCanvas& canvas, TFile& file, const TString& shortName,
                 Int_t first, Int_t last, const char* outputPdf) {
  const char* suffixes[] = {"", "_SBPos", "_SBNeg", "_signalLow", "_signalHigh"};
  const char* labels[] = {"signal", "SBPos", "SBNeg", "signalLow", "signalHigh"};
  for (Int_t region = 0; region < 5; ++region) {
    const TString tail = shortName + suffixes[region] + "_CentBin%d";
    TH2* se = SumMass(file, "true/hKstarMt_" + tail, "qaPairMtSE", first, last);
    TH2* me = SumMass(file, "mix/hKstarMt_Mixed_" + tail, "qaPairMtME", first, last);
    if (!se || !me) {
      std::cerr << "ERROR: could not sum kstar-mT histograms for " << labels[region] << std::endl;
      delete se; delete me; return 1;
    }
    canvas.Clear(); canvas.Divide(2, 1);
    canvas.cd(1);
    se->SetTitle(TString::Format("SE #Lambda-%s %s, cent9 %d-%d;k* (GeV/c);m_{T} (GeV/c^{2})",
                                shortName.Data(), labels[region], first, last));
    se->SetStats(kFALSE); se->Draw("COLZ");
    canvas.cd(2);
    me->SetTitle(TString::Format("ME #Lambda-%s %s, cent9 %d-%d;k* (GeV/c);m_{T} (GeV/c^{2})",
                                shortName.Data(), labels[region], first, last));
    me->SetStats(kFALSE); me->Draw("COLZ");
    canvas.Print(outputPdf);
    delete se; delete me;
  }
  return 0;
}

Int_t CheckHe4RapiditySplit(TFile& file, const TString& channel,
                           const std::string& species, Double_t& threshold) {
  if (species != "he4") return 0;
  const char* regions[] = {"", "_signal", "_leftSB", "_rightSB"};
  const char* sides[] = {"_YLabBelow", "_YLabAbove"};
  Int_t present = 0;
  TH1* reference = dynamic_cast<TH1*>(file.Get("hKstarSE_" + channel));
  TH2* referenceCent = dynamic_cast<TH2*>(file.Get("hKstarSEVsCent_" + channel));
  for (Int_t mixed = 0; mixed < 2; ++mixed) {
    const TString prefix = mixed ? "hKstarME" : "hKstarSE";
    for (Int_t region = 0; region < 4; ++region)
      for (Int_t side = 0; side < 2; ++side)
        for (Int_t dimension = 1; dimension <= 2; ++dimension) {
          const TString name = prefix + (dimension == 2 ? "VsCent_" : "_")
              + channel + regions[region] + sides[side];
          TObject* object = file.Get(name);
          if (!object) continue;
          ++present;
          TH1* histogram = dynamic_cast<TH1*>(object);
          TH2* histogramCent = dynamic_cast<TH2*>(object);
          if (!histogram || histogram->GetDimension() != dimension || !reference
              || !SameAxis(*reference->GetXaxis(), *histogram->GetXaxis())
              || (dimension == 2 && (!histogramCent || !referenceCent
                  || !SameAxis(*referenceCent->GetYaxis(), *histogramCent->GetYaxis())))) {
            std::cerr << "ERROR: incompatible He4 lab-rapidity histogram " << name << std::endl;
            return -1;
          }
        }
  }
  TParameter<double>* cut = dynamic_cast<TParameter<double>*>(
      file.Get("he4RapiditySplitYLabThreshold"));
  if (present == 0 && !cut) return 0;
  if (present != 32 || !cut || !TMath::Finite(cut->GetVal())) {
    std::cerr << "ERROR: incomplete He4 lab-rapidity histogram/threshold contract ("
              << present << "/32)" << std::endl;
    return -1;
  }
  threshold = cut->GetVal();
  return 1;
}

Int_t DrawHe4RapiditySplit(TCanvas& canvas, TFile& file, const TString& channel,
                          Int_t first, Int_t last, Int_t rebin, Double_t threshold,
                          const char* outputPdf) {
  const char* regions[] = {"", "_signal", "_leftSB", "_rightSB"};
  const char* labels[] = {"full mass", "signal", "left sideband", "right sideband"};
  const char* sides[] = {"_YLabBelow", "_YLabAbove"};
  const char* signs[] = {"<", ">"};
  for (Int_t region = 0; region < 4; ++region) {
    canvas.Clear(); canvas.Divide(2, 2);
    TH1* projections[4] = {0, 0, 0, 0};
    for (Int_t side = 0; side < 2; ++side)
      for (Int_t mixed = 0; mixed < 2; ++mixed) {
        const Int_t pad = 2 * side + mixed;
        const TString prefix = mixed ? "hKstarMEVsCent_" : "hKstarSEVsCent_";
        TH2* source = dynamic_cast<TH2*>(file.Get(prefix + channel + regions[region] + sides[side]));
        if (!source) {
          for (Int_t i = 0; i < 4; ++i) delete projections[i];
          return 1;
        }
        source->GetXaxis()->UnZoom();
        projections[pad] = source->ProjectionX(
            TString::Format("qaHe4Y_%d_%d", region, pad),
            source->GetYaxis()->FindBin(first), source->GetYaxis()->FindBin(last), "e");
        projections[pad]->SetDirectory(0);
        projections[pad]->Rebin(rebin);
        projections[pad]->SetStats(kFALSE);
        projections[pad]->SetTitle(TString::Format(
            "%s, ^{4}He Y_{Lab} %s %.3f, %s, cent9 %d-%d;k* (GeV/c);Raw pairs",
            mixed ? "ME" : "SE", signs[side], threshold, labels[region], first, last));
        projections[pad]->SetLineColor(mixed ? kRed + 1 : kBlack);
        canvas.cd(pad + 1); projections[pad]->Draw("E1");
      }
    canvas.Print(outputPdf);
    for (Int_t i = 0; i < 4; ++i) delete projections[i];
  }
  return 0;
}

Int_t CheckEpdQa(TFile& file, TH1* histograms[3]) {
  if (!file.Get("FemtoEpdEventPlaneDefinition")) return 0;
  const char* names[] = {"hEpdPsi1VsCent", "hEpdMixingBinVsCent", "hEpdStatus"};
  const Int_t dimensions[] = {2, 2, 1};
  const Int_t bins[] = {72, 8, 16};
  const Double_t low[] = {-TMath::Pi(), -0.5, -0.5};
  const Double_t high[] = {TMath::Pi(), 7.5, 15.5};
  for (Int_t i = 0; i < 3; ++i) {
    histograms[i] = dynamic_cast<TH1*>(file.Get(names[i]));
    if (!histograms[i] || histograms[i]->GetDimension() != dimensions[i]
        || histograms[i]->GetXaxis()->GetNbins() != bins[i]
        || histograms[i]->GetXaxis()->GetXmin() != low[i]
        || histograms[i]->GetXaxis()->GetXmax() != high[i]
        || (i < 2 && histograms[i]->GetYaxis()->GetNbins() != 9)) {
      std::cerr << "ERROR: incompatible EPD QA histogram " << names[i] << std::endl;
      return -1;
    }
  }
  return 1;
}

Int_t CheckPairRapidity(TFile& file, const TString& channel, const std::string& species) {
  if (species != "he4") return 0;
  const char* regions[] = {"", "_signal", "_leftSB", "_rightSB"};
  Int_t present = 0;
  TParameter<int>* bins = dynamic_cast<TParameter<int>*>(file.Get("he4PairRapidityLabNBins"));
  TParameter<double>* low = dynamic_cast<TParameter<double>*>(file.Get("he4PairRapidityLabMin"));
  TParameter<double>* high = dynamic_cast<TParameter<double>*>(file.Get("he4PairRapidityLabMax"));
  TH2* reference = dynamic_cast<TH2*>(file.Get("hKstarSEVsCent_" + channel));
  for (Int_t mixed = 0; mixed < 2; ++mixed) for (Int_t region = 0; region < 4; ++region) {
    const TString key = TString("hKstar") + (mixed ? "ME" : "SE")
        + "VsPairYVsCent_" + channel + regions[region];
    TObject* object = file.Get(key);
    if (!object) continue;
    ++present;
    TH3D* h = dynamic_cast<TH3D*>(object);
    if (!h || !reference || !bins || !low || !high || bins->GetVal() <= 0
        || !TMath::Finite(low->GetVal()) || !TMath::Finite(high->GetVal())
        || low->GetVal() >= high->GetVal()
        || !SameAxis(*h->GetXaxis(), *reference->GetXaxis())
        || !SameAxis(*h->GetZaxis(), *reference->GetYaxis())
        || h->GetYaxis()->GetNbins() != bins->GetVal()
        || h->GetYaxis()->GetXmin() != low->GetVal()
        || h->GetYaxis()->GetXmax() != high->GetVal()) {
      std::cerr << "ERROR: incompatible pair-rapidity histogram/metadata " << key << std::endl;
      return -1;
    }
  }
  if (!present && !bins && !low && !high && !file.Get("FemtoHe4PairRapiditySplitDefinition")) return 0;
  if (present != 8 || !file.Get("FemtoHe4PairRapiditySplitDefinition")) {
    std::cerr << "ERROR: incomplete pair-rapidity histogram contract " << present << "/8" << std::endl;
    return -1;
  }
  return 1;
}

// Adjacent rapidity bins are merged by specifying a wider inclusive bin range.
// The returned raw k* projection belongs to the caller; flow bins are selectable.
TH1* ProjectPairRapidity(TH3& source, const char* name, Int_t yFirst, Int_t yLast,
                        Int_t centFirst, Int_t centLast) {
  if (yFirst < 0 || yFirst > yLast || yLast > source.GetNbinsY() + 1
      || centFirst < 0 || centLast > 8 || centFirst > centLast) return 0;
  // Native ROOT 5 projection treats (0,0) as the complete axis and can
  // inherit a display zoom. Construct the full X axis without changing source.
  const TAxis& xAxis = *source.GetXaxis();
  TH1D* h = xAxis.GetXbins()->GetSize()
      ? new TH1D(name, source.GetTitle(), xAxis.GetNbins(), xAxis.GetXbins()->GetArray())
      : new TH1D(name, source.GetTitle(), xAxis.GetNbins(), xAxis.GetXmin(), xAxis.GetXmax());
  h->SetDirectory(0); h->Sumw2();
  h->GetXaxis()->SetTitle(xAxis.GetTitle()); h->GetYaxis()->SetTitle("Raw pairs");
  Double_t entries = 0;
  for (Int_t x = 0; x <= source.GetNbinsX() + 1; ++x) {
    Double_t content = 0, variance = 0;
    for (Int_t y = yFirst; y <= yLast; ++y)
      for (Int_t z = source.GetZaxis()->FindBin(centFirst);
           z <= source.GetZaxis()->FindBin(centLast); ++z) {
        const Int_t cell = source.GetBin(x, y, z);
        const Double_t count = source.GetBinContent(cell);
        content += count;
        variance += source.GetSumw2N() ? source.GetSumw2()->At(cell) : TMath::Abs(count);
      }
    h->SetBinContent(x, content); h->SetBinError(x, TMath::Sqrt(variance));
    entries += content;
  }
  h->SetEntries(entries);
  return h;
}

Int_t DrawPairRapidity(TCanvas& canvas, TFile& file, const TString& channel,
                      Int_t first, Int_t last, Int_t rebin, const char* outputPdf) {
  const char* regions[] = {"", "_signal", "_leftSB", "_rightSB"};
  const char* labels[] = {"full mass", "signal", "left sideband", "right sideband"};
  for (Int_t region = 0; region < 4; ++region) for (Int_t mixed = 0; mixed < 2; ++mixed) {
    TH3* source = dynamic_cast<TH3*>(file.Get(TString("hKstar") + (mixed ? "ME" : "SE")
        + "VsPairYVsCent_" + channel + regions[region]));
    if (!source) return 1;
    const Int_t n = source->GetNbinsY();
    // Keep the PDF canvas dimensions fixed and show at most six readable panels.
    for (Int_t start = 1; start <= n; start += 6) {
    canvas.Clear(); canvas.Divide(2, 3);
    std::vector<TH1*> projections(n, 0);
    for (Int_t bin = start; bin <= n && bin < start + 6; ++bin) {
      projections[bin - 1] = ProjectPairRapidity(*source,
          TString::Format("qaPairY_%d_%d_%d", region, mixed, bin), bin, bin, first, last);
      if (!projections[bin - 1]) {
        for (Int_t i = 0; i < n; ++i) delete projections[i];
        return 1;
      }
      TH1* h = projections[bin - 1]; h->Rebin(rebin); h->SetStats(kFALSE);
      h->SetTitle(TString::Format("%s %s, %.2f #leq Y_{pair,Lab} < %.2f;k* (GeV/c);Raw pairs",
          mixed ? "ME" : "SE", labels[region], source->GetYaxis()->GetBinLowEdge(bin),
          source->GetYaxis()->GetBinUpEdge(bin)));
      h->SetLineColor(mixed ? kRed + 1 : kBlack);
      canvas.cd(bin - start + 1); h->Draw("E1");
    }
    canvas.Print(outputPdf);
    for (Int_t i = 0; i < n; ++i) delete projections[i];
    }
    TH1* under = ProjectPairRapidity(*source, "qaPairY_under", 0, 0, first, last);
    TH1* over = ProjectPairRapidity(*source, "qaPairY_over", n + 1, n + 1, first, last);
    std::cout << "[QA pair rapidity] " << (mixed ? "ME " : "SE ") << labels[region]
              << " cent9=" << first << "-" << last
              << " stored Y-underflow=" << (under ? under->Integral(0, under->GetNbinsX()+1) : -1)
              << " Y-overflow=" << (over ? over->Integral(0, over->GetNbinsX()+1) : -1) << std::endl;
    delete under; delete over;
  }
  return 0;
}

void DrawEpdQa(TCanvas& canvas, TH1* histograms[3], const char* outputPdf) {
  canvas.Clear(); canvas.Divide(2, 2);
  for (Int_t i = 0; i < 3; ++i) {
    canvas.cd(i + 1); histograms[i]->SetStats(kFALSE);
    histograms[i]->Draw(i < 2 ? "COLZ" : "HIST");
  }
  canvas.cd(4);
  Message("EP1 and EP8 controls use the same valid EPD events.");
  canvas.Print(outputPdf);
}
}

Int_t checkHistAnaFemtoLambda(const char* rootFile, const char* mainconf, const char* outputPdf) {
  if (!rootFile || !*rootFile || !mainconf || !*mainconf || !outputPdf || !*outputPdf) return 1;
  if (!gSystem->AccessPathName(outputPdf)) {
    std::cerr << "ERROR: refusing existing QA PDF " << outputPdf << std::endl; return 1;
  }
  TString config(mainconf);
  if (!config.BeginsWith("/")) config = TString(gSystem->WorkingDirectory()) + "/" + config;
  std::cout << "[mainconf] (argument) " << config << std::endl;
  if (!ConfigManager::GetInstance().LoadConfig(config)) return 1;
  const FemtoConfig& settings = ConfigManager::GetInstance().GetFemtoConfig();
  if (settings.channels.empty() || settings.channels[0].partA != "lambda") return 1;
  const FemtoConfig::ChannelDef& base = settings.channels[0];
  const std::string species = base.partB;
  TString shortName;
  if (species == "deuteron") shortName = "d";
  else if (species == "triton") shortName = "t";
  else if (species == "he3") shortName = "3He";
  else if (species == "he4") shortName = "4He";
  else return 1;
  TFile file(rootFile, "READ");
  if (file.IsZombie()) return 1;
  TNamed* status = dynamic_cast<TNamed*>(file.Get("FemtoLambdaRunStatus"));
  if (status && TString(status->GetTitle()) != "completed") {
    std::cerr << "ERROR: input run was not completed" << std::endl; return 1;
  }
  const Int_t first = settings.cfCent9Min, last = settings.cfCent9Max;
  if (first < 0 || last > 8 || first > last || settings.cfRebinFactor < 1) return 1;
  TH1* inclusive = dynamic_cast<TH1*>(file.Get("hLambda_InvMass"));
  TH2* massSE = FemtoLambdaQa::SumMass(file, "true/hKstarMass_" + shortName + "_CentBin%d", "qaMassSE", first, last);
  TH2* massME = FemtoLambdaQa::SumMass(file, "mix/hKstarMass_Mixed_" + shortName + "_CentBin%d", "qaMassME", first, last);
  if (!inclusive || !massSE || !massME) {
    std::cerr << "ERROR: missing Lambda mass or full-mass legacy histogram contract" << std::endl;
    delete massSE; delete massME; return 1;
  }
  TH2* acceptance[4] = {0, 0, 0, 0};
  const Int_t extensions = FemtoLambdaQa::CheckAcceptanceAndSplit(
      file, base.name.c_str(), species.c_str(), acceptance);
  if (extensions < 0) { delete massSE; delete massME; return 1; }
  const Int_t pairMt = FemtoLambdaQa::CheckPairMt(file, shortName);
  if (pairMt < 0) { delete massSE; delete massME; return 1; }
  Double_t yLabThreshold = 0;
  const Int_t yLabSplit = FemtoLambdaQa::CheckHe4RapiditySplit(
      file, base.name.c_str(), species, yLabThreshold);
  if (yLabSplit < 0) { delete massSE; delete massME; return 1; }
  const Int_t pairRapidity = FemtoLambdaQa::CheckPairRapidity(file, base.name.c_str(), species);
  if (pairRapidity < 0) { delete massSE; delete massME; return 1; }
  TH1* epdHistograms[3] = {0, 0, 0};
  const Int_t epdQa = FemtoLambdaQa::CheckEpdQa(file, epdHistograms);
  if (epdQa < 0) { delete massSE; delete massME; return 1; }
  const TString outputDirectory = gSystem->DirName(outputPdf);
  if (gSystem->mkdir(outputDirectory.Data(), kTRUE) != 0
      && gSystem->AccessPathName(outputDirectory.Data())) {
    std::cerr << "ERROR: cannot create QA output directory " << outputDirectory << std::endl;
    delete massSE; delete massME; return 1;
  }
  TCanvas canvas("lambdaFemtoQa", "Lambda-nucleus QA", 1200, 900);
  canvas.Print(TString(outputPdf) + "[");
  canvas.Divide(2, 2);
  canvas.cd(1); inclusive->SetTitle("Inclusive selected KF Lambda;M(p#pi) (GeV/c^{2});Candidates"); inclusive->Draw("E");
  canvas.cd(2); massSE->SetTitle("SE pair Lambda mass vs k*;k* (GeV/c);M(p#pi) (GeV/c^{2})"); massSE->Draw("COLZ");
  canvas.cd(3); massME->SetTitle("ME pair Lambda mass vs k*;k* (GeV/c);M(p#pi) (GeV/c^{2})"); massME->Draw("COLZ");
  canvas.cd(4);
  TH1* pairMass = massSE->ProjectionY("qaPairLambdaMass", 0, massSE->GetNbinsX()+1);
  pairMass->SetDirectory(0);
  pairMass->SetTitle("Pair-associated Lambda mass (SE, all k* including flow);M(p#pi) (GeV/c^{2});Pairs");
  pairMass->Draw("E");
  canvas.Print(outputPdf);
  delete pairMass;

  Int_t result = 0;
  const char* suffixes[] = {"signal", "leftSB", "rightSB", "signalLow", "signalHigh"};
  const Int_t regions = extensions ? 5 : 3;
  for (Int_t region = 0; region < regions; ++region) {
    // Preserve the original four pages first, then add acceptance and the two split CFs.
    if (region == 3) FemtoLambdaQa::DrawAcceptance(canvas, acceptance, outputPdf);
    canvas.Clear(); canvas.Divide(2, 1);
    const TString channel = TString(base.name.c_str()) + "_" + suffixes[region];
    TH2* seCent = dynamic_cast<TH2*>(file.Get("hKstarSEVsCent_" + channel));
    TH2* meCent = dynamic_cast<TH2*>(file.Get("hKstarMEVsCent_" + channel));
    if (!seCent || !meCent) {
      std::cerr << "ERROR: missing SE/ME centrality histogram for " << channel << std::endl;
      result = 1; break;
    }
    TH1* se = seCent->ProjectionX("qaSE", seCent->GetYaxis()->FindBin(first), seCent->GetYaxis()->FindBin(last));
    TH1* me = meCent->ProjectionX("qaME", meCent->GetYaxis()->FindBin(first), meCent->GetYaxis()->FindBin(last));
    se->SetDirectory(0); me->SetDirectory(0);
    se->Rebin(settings.cfRebinFactor); me->Rebin(settings.cfRebinFactor);
    const Int_t normFirst = se->GetXaxis()->FindFixBin(base.normQMin);
    const Int_t normLast = se->GetXaxis()->FindFixBin(base.normQMax - 1e-9);
    const Double_t seNorm = se->Integral(normFirst, normLast), meNorm = me->Integral(normFirst, normLast);
    canvas.cd(1);
    se->SetTitle(TString::Format("%s, cent9 %d-%d;k* (GeV/c);Pairs", channel.Data(), first, last));
    se->SetLineColor(kBlack); me->SetLineColor(kRed+1);
    se->SetStats(kFALSE); me->SetStats(kFALSE);
    if (seNorm > 0 && meNorm > 0) me->Scale(seNorm / meNorm);
    se->SetMaximum(TMath::Max(1.0, 1.15 * TMath::Max(se->GetMaximum(), me->GetMaximum())));
    se->Draw("E"); me->Draw("HIST SAME");
    TLegend legend(0.53, 0.73, 0.89, 0.89); legend.AddEntry(se, "SE", "l");
    legend.AddEntry(me, seNorm > 0 && meNorm > 0 ? "ME (norm region scaled)" : "ME (unscaled)", "l"); legend.Draw();
    canvas.cd(2);
    TH1* cf = 0;
    TGraphErrors validCf;
    validCf.SetMarkerStyle(20);
    if (seNorm > 0 && meNorm > 0) {
      cf = static_cast<TH1*>(se->Clone("qaRawCF")); cf->SetDirectory(0); cf->Divide(me);
      cf->SetTitle("Raw CF (zero-ME bins omitted); k* (GeV/c);C(k*)");
      // Clone inherits SE display bounds; reset before querying ratio data maximum.
      cf->SetMaximum(); cf->SetStats(kFALSE);
      cf->SetMinimum(0); cf->SetMaximum(TMath::Max(1.5, 1.15*cf->GetMaximum()));
      cf->Draw("AXIS");
      for (Int_t bin = 1; bin <= cf->GetNbinsX(); ++bin) {
        if (me->GetBinContent(bin) <= 0) continue;
        const Int_t point = validCf.GetN();
        validCf.SetPoint(point, cf->GetBinCenter(bin), cf->GetBinContent(bin));
        validCf.SetPointError(point, 0, cf->GetBinError(bin));
      }
      validCf.Draw("P SAME");
      std::cout << "[QA] " << channel << " normalization SE=" << seNorm << " ME=" << meNorm << std::endl;
    } else {
      FemtoLambdaQa::Message("Insufficient normalization counts; CF is not defined.");
      std::cout << "[QA] " << channel << " CF undefined: zero normalization counts" << std::endl;
    }
    canvas.Print(outputPdf);
    delete cf; delete se; delete me;
  }
  if (result == 0 && pairMt)
    result = FemtoLambdaQa::DrawPairMt(canvas, file, shortName, first, last, outputPdf);
  if (result == 0 && yLabSplit)
    result = FemtoLambdaQa::DrawHe4RapiditySplit(canvas, file, base.name.c_str(), first, last,
                                               settings.cfRebinFactor, yLabThreshold, outputPdf);
  if (result == 0 && epdQa)
    FemtoLambdaQa::DrawEpdQa(canvas, epdHistograms, outputPdf);
  if (result == 0 && pairRapidity)
    result = FemtoLambdaQa::DrawPairRapidity(canvas, file, base.name.c_str(), first, last,
                                           settings.cfRebinFactor, outputPdf);
  canvas.Print(TString(outputPdf) + "]");
  delete massSE; delete massME;
  FileStat_t pdfInfo;
  if (gSystem->GetPathInfo(outputPdf, pdfInfo) != 0 || pdfInfo.fSize <= 0) {
    std::cerr << "ERROR: QA PDF was not written" << std::endl; return 1;
  }
  return result;
}
