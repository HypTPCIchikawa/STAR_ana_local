// ROOT 5 CINT needs only the public entry point. ACLiC compiles the helpers.
void make_cf_plots(const char* binsTsv, const char* metaTsv, const char* panelsTsv,
                   const char* summaryTsv, const char* outDir, const char* outRoot);

#ifndef __CINT__

#include "TArrayD.h"
#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TObject.h"
#include "TPad.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ep_rapidity_cf_plots {

bool rootErrorSeen = false;
std::string rootErrorText;
void RootErrorHandler(Int_t level, Bool_t, const char* location, const char* message) {
  if (level >= kError) {
    rootErrorSeen = true;
    rootErrorText = std::string(location ? location : "ROOT") + ": " + (message ? message : "");
  }
  DefaultErrorHandler(level, kFALSE, location, message);
}
struct RootErrorGuard {
  ErrorHandlerFunc_t previous;
  RootErrorGuard() : previous(SetErrorHandler(RootErrorHandler)) { rootErrorSeen = false; rootErrorText.clear(); }
  ~RootErrorGuard() { SetErrorHandler(previous); }
};
void Require(bool good, const std::string& message) { if (!good) throw std::runtime_error(message); }
void CheckRoot() { Require(!rootErrorSeen, std::string("ROOT error: ") + rootErrorText); }
bool Finite(double value) { return TMath::Finite(value) != 0; }
double Number(const std::string& text) {
  std::istringstream stream(text); double value = 0.;
  Require(bool(stream >> value) && Finite(value), std::string("Invalid/nonfinite number: ") + text);
  stream >> std::ws; Require(stream.eof(), std::string("Trailing numeric text: ") + text); return value;
}
int Integer(const std::string& text) {
  std::istringstream stream(text); int value = 0;
  Require(bool(stream >> value), std::string("Invalid integer: ") + text);
  stream >> std::ws; Require(stream.eof(), std::string("Trailing integer text: ") + text); return value;
}
struct Optional {
  bool valid; double value;
  Optional() : valid(false), value(0.) {}
  explicit Optional(const std::string& text) : valid(text != "NA"), value(0.) { if (valid) value = Number(text); }
};
bool SafeName(const std::string& text) {
  if (text.empty() || text == "." || text == "..") return false;
  for (size_t i = 0; i < text.size(); ++i) {
    const char c = text[i];
    if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
          (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.')) return false;
  }
  return true;
}
std::string Key(const std::string& dataset, const std::string& cent, const std::string& selection) {
  Require((dataset == "EP1" || dataset == "EP8") && SafeName(cent) && SafeName(selection), "Invalid dataset/centrality/selection identifier");
  return dataset + "\t" + cent + "\t" + selection;
}
typedef std::vector<std::string> Fields;
typedef std::vector<Fields> Table;
void ReadTable(const char* path, size_t columns, Table& table) {
  Require(path && *path, "An input TSV path is empty");
  std::ifstream stream(path); Require(bool(stream), std::string("Cannot read TSV: ") + path);
  std::string line; unsigned int lineNumber = 0;
  while (std::getline(stream, line)) {
    ++lineNumber;
    if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
    if (line.empty()) continue;
    Fields fields; size_t start = 0;
    for (;;) {
      const size_t end = line.find('\t', start);
      fields.push_back(line.substr(start, end == std::string::npos ? end : end - start));
      if (end == std::string::npos) break;
      start = end + 1;
    }
    std::ostringstream context; context << path << ':' << lineNumber;
    Require(fields.size() == columns, context.str() + ": wrong TSV column count");
    table.push_back(fields);
  }
  Require(!stream.bad() && !table.empty(), std::string("Empty/failed TSV: ") + path);
}

struct Bin {
  int index;
  double low, high, se, seVariance, me, meVariance;
  Optional cf, cfError;
};
struct Meta {
  bool present; Optional scale;
  double normLow, normHigh, displayLow, displayHigh, highLow, highHigh;
  Meta() : present(false), normLow(0), normHigh(0), displayLow(0), displayHigh(0), highLow(0), highHigh(0) {}
};
struct Selection {
  std::string dataset, cent, name;
  std::vector<Bin> bins;
  Meta meta;
};
typedef std::map<std::string, Selection> Selections;
bool BinLess(const Bin& a, const Bin& b) { return a.low < b.low; }
void ReadSelections(const char* binsPath, const char* metaPath, Selections& selections) {
  Table metaRows, binRows; ReadTable(metaPath, 10, metaRows); ReadTable(binsPath, 12, binRows);
  for (size_t i = 0; i < metaRows.size(); ++i) {
    const Fields& f = metaRows[i]; Selection& s = selections[Key(f[0], f[1], f[2])];
    Require(!s.meta.present, "Duplicate metadata selection");
    s.dataset = f[0]; s.cent = f[1]; s.name = f[2]; s.meta.present = true; s.meta.scale = Optional(f[3]);
    s.meta.normLow = Number(f[4]); s.meta.normHigh = Number(f[5]);
    s.meta.displayLow = Number(f[6]); s.meta.displayHigh = Number(f[7]);
    s.meta.highLow = Number(f[8]); s.meta.highHigh = Number(f[9]);
    Require(s.meta.normLow < s.meta.normHigh && s.meta.displayLow < s.meta.displayHigh && s.meta.highLow < s.meta.highHigh,
            "Reversed metadata interval");
    Require(!s.meta.scale.valid || s.meta.scale.value > 0., "Invalid ME scale");
  }
  for (size_t i = 0; i < binRows.size(); ++i) {
    const Fields& f = binRows[i]; const std::string key = Key(f[0], f[1], f[2]);
    Require(selections.count(key) != 0, "Bin selection is absent from metadata");
    Bin bin; bin.index = Integer(f[3]); bin.low = Number(f[4]); bin.high = Number(f[5]);
    bin.se = Number(f[6]); bin.seVariance = Number(f[7]); bin.me = Number(f[8]); bin.meVariance = Number(f[9]);
    bin.cf = Optional(f[10]); bin.cfError = Optional(f[11]);
    Require(bin.low < bin.high && bin.seVariance >= 0. && bin.meVariance >= 0., "Invalid bin extent/variance");
    Require(bin.cf.valid == bin.cfError.valid && (!bin.cfError.valid || bin.cfError.value >= 0.), "Invalid CF/error availability");
    Require(!bin.cf.valid || selections[key].meta.scale.valid, "CF is supplied without a valid normalization");
    selections[key].bins.push_back(bin);
  }
  for (Selections::iterator it = selections.begin(); it != selections.end(); ++it) {
    Selection& s = it->second; Require(!s.bins.empty(), "Metadata selection has no bin rows");
    std::sort(s.bins.begin(), s.bins.end(), BinLess); std::set<int> indices;
    for (size_t i = 0; i < s.bins.size(); ++i) {
      Require(indices.insert(s.bins[i].index).second, "Duplicate bin index");
      if (i) Require(s.bins[i - 1].high == s.bins[i].low, "Histogram bin edges are not exactly contiguous");
    }
  }
}

struct Panel { std::string cent, selection, title; };
void ReadPanels(const char* path, std::vector<Panel>& panels) {
  Table table; ReadTable(path, 3, table); std::set<std::string> unique;
  for (size_t i = 0; i < table.size(); ++i) {
    Panel p; p.cent = table[i][0]; p.selection = table[i][1]; p.title = table[i][2];
    Require(SafeName(p.cent) && SafeName(p.selection) && !p.title.empty(), "Invalid panel row");
    Require(unique.insert(p.cent + "\t" + p.selection).second, "Duplicate panel"); panels.push_back(p);
  }
}
Selection& Get(Selections& selections, const std::string& dataset, const Panel& panel) {
  const std::string key = Key(dataset, panel.cent, panel.selection);
  Require(selections.count(key) != 0, "Panel selection is absent from inputs"); return selections.find(key)->second;
}
void SamePanel(const Selection& a, const Selection& b) {
  Require(a.bins.size() == b.bins.size(), "EP1/EP8 panel bin counts differ");
  for (size_t i = 0; i < a.bins.size(); ++i)
    Require(a.bins[i].low == b.bins[i].low && a.bins[i].high == b.bins[i].high, "EP1/EP8 panel bin edges differ");
  Require(a.meta.normLow == b.meta.normLow && a.meta.normHigh == b.meta.normHigh &&
          a.meta.displayLow == b.meta.displayLow && a.meta.displayHigh == b.meta.displayHigh &&
          a.meta.highLow == b.meta.highLow && a.meta.highHigh == b.meta.highHigh, "EP1/EP8 display/normalization intervals differ");
}

int Color(const std::string& dataset) { return TColor::GetColor(dataset == "EP1" ? "#244C75" : "#C05A24"); }
int Marker(const std::string& dataset) { return dataset == "EP1" ? 20 : 21; }
void GraphStyle(TGraphErrors& graph, const std::string& dataset) {
  graph.SetLineColor(Color(dataset)); graph.SetMarkerColor(Color(dataset));
  graph.SetMarkerStyle(Marker(dataset)); graph.SetMarkerSize(0.85); graph.SetLineWidth(2);
}
void CFGraph(const Selection& s, TGraphErrors& graph) {
  graph.SetName("CF_graph"); graph.SetTitle(";k^{*} (GeV/c);C(k^{*})");
  for (size_t i = 0; i < s.bins.size(); ++i) {
    const Bin& bin = s.bins[i]; if (!bin.cf.valid) continue;
    const int point = graph.GetN();
    graph.SetPoint(point, (bin.low + bin.high) / 2., bin.cf.value);
    graph.SetPointError(point, (bin.high - bin.low) / 2., bin.cfError.value);
  }
  GraphStyle(graph, s.dataset);
}
void Histogram(const Selection& s, bool same, const std::string& name, TH1D*& histogram) {
  std::vector<double> edges(s.bins.size() + 1);
  for (size_t i = 0; i < s.bins.size(); ++i) edges[i] = s.bins[i].low;
  edges[s.bins.size()] = s.bins.back().high;
  histogram = new TH1D(name.c_str(), ";k^{*} (GeV/c);Pair counts / bin", int(s.bins.size()), &edges[0]);
  histogram->SetDirectory(0); histogram->Sumw2(); double entries = 0.;
  for (size_t i = 0; i < s.bins.size(); ++i) {
    const Bin& b = s.bins[i]; const double value = same ? b.se : b.me;
    histogram->SetBinContent(int(i) + 1, value);
    histogram->GetSumw2()->SetAt(same ? b.seVariance : b.meVariance, int(i) + 1); entries += value;
  }
  Require(Finite(entries), "Nonfinite histogram entries"); histogram->SetEntries(entries);
}
TDirectory* Directory(TDirectory* parent, const std::string& name) {
  TDirectory* directory = parent->GetDirectory(name.c_str()); if (!directory) directory = parent->mkdir(name.c_str());
  CheckRoot(); Require(directory != 0, "Cannot create ROOT directory"); return directory;
}
void WriteDouble(const char* name, double value) {
  TParameter<double> parameter(name, value); Require(parameter.Write() > 0, "Cannot write numeric metadata");
}
void WriteSelection(TFile* output, const Selection& s) {
  TDirectory* directory = Directory(Directory(Directory(output, s.dataset), s.cent), s.name); directory->cd();
  TH1D *se = 0, *me = 0; Histogram(s, true, "SE", se); Histogram(s, false, "ME", me);
  TGraphErrors cf; CFGraph(s, cf);
  Require(se->Write() > 0 && me->Write() > 0 && cf.Write("CF_graph") > 0, "Cannot write raw histograms/CF graph");
  if (s.meta.scale.valid && cf.GetN() > 0) {
    TH1D scaled(*me); scaled.SetName("ME_scaled"); scaled.SetDirectory(0); scaled.Scale(s.meta.scale.value);
    Require(scaled.Write() > 0, "Cannot write scaled ME"); WriteDouble("normA", s.meta.scale.value);
  }
  WriteDouble("normLow", s.meta.normLow); WriteDouble("normHigh", s.meta.normHigh);
  WriteDouble("displayLow", s.meta.displayLow); WriteDouble("displayHigh", s.meta.displayHigh);
  WriteDouble("highLow", s.meta.highLow); WriteDouble("highHigh", s.meta.highHigh);
  TParameter<int> valid("nValidCFBins", cf.GetN()); Require(valid.Write() > 0, "Cannot write valid-bin count");
  delete se; delete me; CheckRoot();
}

struct Range {
  double low, high;
  Range() : low(1.), high(1.) {}
  void Add(double value, double error) { low = std::min(low, value - error); high = std::max(high, value + error); }
  void Pad() {
    double span = high - low; if (!(span > 0.)) span = 0.2;
    low -= 0.14 * span; high += 0.22 * span;
    Require(Finite(low) && Finite(high) && low < high, "Invalid plot y range");
  }
};
void ExtendCF(const Selection& s, double low, double high, Range& range) {
  for (size_t i = 0; i < s.bins.size(); ++i) {
    const Bin& b = s.bins[i];
    if (b.cf.valid && b.high > low && b.low < high) range.Add(b.cf.value, b.cfError.value);
  }
}
void PadStyle(TPad* pad) {
  pad->SetLeftMargin(0.13); pad->SetRightMargin(0.04); pad->SetTopMargin(0.11); pad->SetBottomMargin(0.14);
  pad->SetTicks(1, 1); pad->SetGridy(kFALSE);
}
void Text(double x, double y, const std::string& text, double size = 0.035) {
  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(size); latex.DrawLatex(x, y, text.c_str());
}
std::string NumberText(double value) { std::ostringstream text; text << std::setprecision(4) << value; return text.str(); }
void Save(TCanvas& canvas, const std::string& path) {
  canvas.Print(path.c_str()); CheckRoot(); std::ifstream file(path.c_str(), std::ios::binary | std::ios::ate);
  Require(bool(file) && file.tellg() > 0, std::string("Missing/empty PDF: ") + path);
}
void CFPanel(const Panel& panel, const Selection& ep1, const Selection& ep8, bool highOnly, const std::string& path) {
  const double xLow = highOnly ? ep1.meta.highLow : ep1.meta.displayLow;
  const double xHigh = highOnly ? ep1.meta.highHigh : ep1.meta.displayHigh;
  Range range; ExtendCF(ep1, xLow, xHigh, range); ExtendCF(ep8, xLow, xHigh, range); range.Pad();
  TCanvas canvas("cf_canvas", "CF comparison", 1000, 560); PadStyle(&canvas);
  TH1D frame("cf_frame", ";k^{*} (GeV/c);C(k^{*})", 1, xLow, xHigh); frame.SetDirectory(0);
  frame.SetStats(kFALSE); frame.SetMinimum(range.low); frame.SetMaximum(range.high);
  frame.GetXaxis()->SetTitleSize(0.045); frame.GetYaxis()->SetTitleSize(0.045); frame.GetYaxis()->SetTitleOffset(1.25);
  frame.Draw("AXIS");
  const double normLow = std::max(xLow, ep1.meta.normLow), normHigh = std::min(xHigh, ep1.meta.normHigh);
  TBox band(normLow, range.low, normHigh, range.high); band.SetFillColor(TColor::GetColor("#E5E7E9")); band.SetLineColor(0);
  if (normLow < normHigh) band.Draw("same");
  TLine unity(xLow, 1., xHigh, 1.); unity.SetLineStyle(2); unity.SetLineColor(kGray + 2); unity.Draw("same");
  TGraphErrors g1, g8; CFGraph(ep1, g1); CFGraph(ep8, g8);
  if (g1.GetN()) g1.Draw("P same"); if (g8.GetN()) g8.Draw("P same");
  TLegend legend(0.67, 0.74, 0.94, 0.87); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42);
  legend.AddEntry(&g1, "EP1", "pe"); legend.AddEntry(&g8, "EP8", "pe"); legend.Draw();
  Text(0.13, 0.94, panel.title + (highOnly ? " : high k^{*}" : ""), 0.039);
  Text(0.15, 0.83, "Normalization: " + NumberText(ep1.meta.normLow) + " < k^{*} < " + NumberText(ep1.meta.normHigh) + " GeV/c", 0.029);
  if (!g1.GetN() || !g8.GetN()) Text(0.15, 0.77, "Unavailable CF points are omitted", 0.028);
  canvas.RedrawAxis(); Save(canvas, path); canvas.Clear();
}
void SpectraPanel(const Panel& panel, const Selection& ep1, const Selection& ep8, const std::string& path) {
  TCanvas canvas("spectra_canvas", "SE and normalized ME", 1000, 850); canvas.Divide(1, 2, 0., 0.005);
  const Selection* data[] = {&ep1, &ep8};
  TH1D *same[2] = {0, 0}, *mixed[2] = {0, 0};
  for (int d = 0; d < 2; ++d) {
    canvas.cd(d + 1); PadStyle(static_cast<TPad*>(canvas.GetPad(d + 1)));
    const Selection& s = *data[d]; Histogram(s, true, "spectra_SE_" + s.dataset, same[d]);
    Histogram(s, false, "spectra_ME_" + s.dataset, mixed[d]);
    const bool scaled = s.meta.scale.valid; if (scaled) mixed[d]->Scale(s.meta.scale.value);
    double maximum = 0.;
    for (size_t i = 0; i < s.bins.size(); ++i) if (s.bins[i].high > s.meta.displayLow && s.bins[i].low < s.meta.displayHigh) {
      maximum = std::max(maximum, same[d]->GetBinContent(int(i) + 1) + same[d]->GetBinError(int(i) + 1));
      if (scaled) maximum = std::max(maximum, mixed[d]->GetBinContent(int(i) + 1));
    }
    if (!(maximum > 0.)) maximum = 1.;
    same[d]->SetStats(kFALSE); same[d]->SetMinimum(0.); same[d]->SetMaximum(maximum * 1.25);
    same[d]->GetXaxis()->SetRangeUser(s.meta.displayLow, s.meta.displayHigh);
    same[d]->SetLineColor(kBlack); same[d]->SetMarkerColor(kBlack); same[d]->SetMarkerStyle(20); same[d]->SetMarkerSize(0.65);
    same[d]->Draw("E1");
    mixed[d]->SetLineColor(TColor::GetColor("#B32628")); mixed[d]->SetLineWidth(2);
    if (scaled) mixed[d]->Draw("HIST same");
    TLegend legend(0.65, 0.72, 0.95, 0.88); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42);
    legend.AddEntry(same[d], "Same event", "pe");
    if (scaled) legend.AddEntry(mixed[d], "Mixed event (scaled)", "l"); legend.DrawClone();
    Text(0.13, 0.94, panel.title + " : " + s.dataset, 0.038);
    if (scaled) Text(0.15, 0.82, "ME scale = " + NumberText(s.meta.scale.value), 0.03);
    else Text(0.15, 0.82, "Normalization unavailable; SE shown", 0.03);
    gPad->RedrawAxis();
  }
  Save(canvas, path); canvas.Clear(); for (int d = 0; d < 2; ++d) { delete same[d]; delete mixed[d]; }
}

struct SummaryPoint {
  std::string dataset, cent, family, label;
  double x, xError; Optional tail, tailError, shape, shapeError;
};
typedef std::map<std::string, std::vector<SummaryPoint> > Summaries;
bool PointLess(const SummaryPoint& a, const SummaryPoint& b) { return a.x < b.x; }
void ReadSummaries(const char* path, Summaries& summaries) {
  Table table; ReadTable(path, 10, table); std::set<std::string> unique;
  for (size_t i = 0; i < table.size(); ++i) {
    const Fields& f = table[i]; SummaryPoint p;
    p.dataset = f[0]; p.cent = f[1]; p.family = f[2];
    Require((p.dataset == "EP1" || p.dataset == "EP8") && SafeName(p.cent) &&
            (p.family == "pair10" || p.family == "pair5" || p.family == "he4"), "Invalid summary identifier");
    p.x = Number(f[3]); p.xError = Number(f[4]); p.label = f[5];
    p.tail = Optional(f[6]); p.tailError = Optional(f[7]); p.shape = Optional(f[8]); p.shapeError = Optional(f[9]);
    Require(p.xError >= 0. && p.tail.valid == p.tailError.valid && p.shape.valid == p.shapeError.valid &&
            (!p.tailError.valid || p.tailError.value >= 0.) && (!p.shapeError.valid || p.shapeError.value >= 0.), "Invalid summary errors");
    Require(unique.insert(p.dataset + "\t" + p.cent + "\t" + p.family + "\t" + NumberText(p.x)).second, "Duplicate summary point");
    summaries[p.cent + "\t" + p.family].push_back(p);
  }
  for (Summaries::iterator it = summaries.begin(); it != summaries.end(); ++it)
    std::sort(it->second.begin(), it->second.end(), PointLess);
}
void SummaryPanel(const std::vector<SummaryPoint>& points, bool shape, const std::string& path) {
  Require(!points.empty(), "Empty summary panel");
  const bool categorical = points[0].family == "he4"; const std::string cent = points[0].cent, family = points[0].family;
  double xLow = points[0].x - points[0].xError, xHigh = points[0].x + points[0].xError;
  Range range; std::map<double, std::string> labels; TGraphErrors ep1, ep8;
  GraphStyle(ep1, "EP1"); GraphStyle(ep8, "EP8");
  for (size_t i = 0; i < points.size(); ++i) {
    const SummaryPoint& p = points[i];
    xLow = std::min(xLow, p.x - p.xError); xHigh = std::max(xHigh, p.x + p.xError);
    if (labels.count(p.x)) Require(labels[p.x] == p.label, "EP1/EP8 category labels differ"); else labels[p.x] = p.label;
    const Optional& y = shape ? p.shape : p.tail; const Optional& error = shape ? p.shapeError : p.tailError;
    if (!y.valid) continue;
    TGraphErrors& graph = p.dataset == "EP1" ? ep1 : ep8; const int n = graph.GetN();
    graph.SetPoint(n, p.x, y.value); graph.SetPointError(n, p.xError, error.value); range.Add(y.value, error.value);
  }
  if (categorical) { xLow = labels.begin()->first - 0.5; xHigh = labels.rbegin()->first + 0.5; }
  else if (xLow == xHigh) { xLow -= 0.5; xHigh += 0.5; }
  else { const double padding = (xHigh - xLow) * 0.025; xLow -= padding; xHigh += padding; }
  range.Pad(); TCanvas canvas("summary_canvas", "High-kstar summary", 1000, 560); PadStyle(&canvas);
  if (categorical) canvas.SetBottomMargin(0.20);
  const std::string title = std::string(";") + (categorical ? "^{4}He y_{lab} selection" : "y_{pair}^{lab}") +
      ";" + (shape ? "High-k^{*} shape ratio" : "High-k^{*} normalized CF");
  TH1D frame("summary_frame", title.c_str(), categorical ? int(labels.size()) : 1, xLow, xHigh);
  frame.SetDirectory(0); frame.SetStats(kFALSE); frame.SetMinimum(range.low); frame.SetMaximum(range.high);
  if (categorical) for (std::map<double, std::string>::const_iterator it = labels.begin(); it != labels.end(); ++it)
    frame.GetXaxis()->SetBinLabel(frame.GetXaxis()->FindBin(it->first), it->second.c_str());
  frame.GetYaxis()->SetTitleOffset(1.3); frame.Draw("AXIS");
  TLine unity(xLow, 1., xHigh, 1.); unity.SetLineStyle(2); unity.SetLineColor(kGray + 2); unity.Draw("same");
  if (ep1.GetN()) ep1.Draw("P same"); if (ep8.GetN()) ep8.Draw("P same");
  TLegend legend(0.70, 0.75, 0.94, 0.88); legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.AddEntry(&ep1, "EP1", "pe"); legend.AddEntry(&ep8, "EP8", "pe"); legend.Draw();
  Text(0.13, 0.94, "Centrality " + cent + " : " + family + (shape ? " : shape" : " : tail"), 0.039);
  if (!ep1.GetN() || !ep8.GetN()) Text(0.15, 0.83, "Unavailable values are omitted", 0.028);
  canvas.RedrawAxis(); Save(canvas, path); canvas.Clear();
}

struct RootOutput {
  TFile* file;
  explicit RootOutput(const char* path) : file(TFile::Open(path, "CREATE")) {
    if (rootErrorSeen || !file || file->IsZombie() || !file->IsWritable()) {
      if (file) { file->Close(); delete file; file = 0; } CheckRoot(); throw std::runtime_error("Cannot create plot ROOT");
    }
  }
  ~RootOutput() { if (file) { file->Close(); delete file; } }
  void Close() { if (file) { file->Flush(); file->Close(); delete file; file = 0; } CheckRoot(); }
};
void Run(const char* binsTsv, const char* metaTsv, const char* panelsTsv, const char* summaryTsv,
         const char* outDir, const char* outRoot) {
  RootErrorGuard errors; const Bool_t previousAddDirectory = TH1::AddDirectoryStatus(); TH1::AddDirectory(kFALSE);
  struct AddDirectoryGuard { Bool_t previous; explicit AddDirectoryGuard(Bool_t p) : previous(p) {} ~AddDirectoryGuard() { TH1::AddDirectory(previous); } } guard(previousAddDirectory);
  gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(3);
  Require(outDir && *outDir && outRoot && *outRoot, "Output paths are required");
  Require(gSystem->AccessPathName(outRoot), "Refusing to overwrite existing output ROOT");
  Selections selections; std::vector<Panel> panels; Summaries summaries;
  ReadSelections(binsTsv, metaTsv, selections); ReadPanels(panelsTsv, panels); ReadSummaries(summaryTsv, summaries);
  std::set<std::string> planned;
  for (size_t i = 0; i < panels.size(); ++i) {
    const Panel& p = panels[i]; SamePanel(Get(selections, "EP1", p), Get(selections, "EP8", p));
    const std::string stem = std::string(outDir) + "/" + p.cent + "__" + p.selection;
    Require(planned.insert(stem + "__cf.pdf").second && planned.insert(stem + "__high.pdf").second &&
            planned.insert(stem + "__spectra.pdf").second, "Duplicate plot target");
  }
  for (Summaries::const_iterator it = summaries.begin(); it != summaries.end(); ++it) {
    const SummaryPoint& p = it->second[0]; const std::string stem = std::string(outDir) + "/" + p.cent + "__" + p.family;
    Require(planned.insert(stem + "__tail_summary.pdf").second && planned.insert(stem + "__shape_summary.pdf").second, "Duplicate summary target");
  }
  for (std::set<std::string>::const_iterator it = planned.begin(); it != planned.end(); ++it)
    Require(gSystem->AccessPathName(it->c_str()), std::string("Refusing to overwrite plot: ") + *it);
  if (gSystem->AccessPathName(outDir)) Require(gSystem->mkdir(outDir, kTRUE) == 0, "Cannot create plot directory");
  RootOutput output(outRoot);
  for (Selections::const_iterator it = selections.begin(); it != selections.end(); ++it) WriteSelection(output.file, it->second);
  for (size_t i = 0; i < panels.size(); ++i) {
    const Panel& p = panels[i]; const Selection& ep1 = Get(selections, "EP1", p); const Selection& ep8 = Get(selections, "EP8", p);
    const std::string stem = std::string(outDir) + "/" + p.cent + "__" + p.selection;
    CFPanel(p, ep1, ep8, false, stem + "__cf.pdf"); CFPanel(p, ep1, ep8, true, stem + "__high.pdf");
    SpectraPanel(p, ep1, ep8, stem + "__spectra.pdf");
  }
  for (Summaries::const_iterator it = summaries.begin(); it != summaries.end(); ++it) {
    const SummaryPoint& p = it->second[0]; const std::string stem = std::string(outDir) + "/" + p.cent + "__" + p.family;
    SummaryPanel(it->second, false, stem + "__tail_summary.pdf"); SummaryPanel(it->second, true, stem + "__shape_summary.pdf");
  }
  output.Close(); CheckRoot();
  std::cout << "MAKE_CF_PLOTS_OK panels=" << panels.size() << " summary_groups=" << summaries.size() << std::endl;
}

} // namespace ep_rapidity_cf_plots

void make_cf_plots(const char* binsTsv, const char* metaTsv, const char* panelsTsv,
                   const char* summaryTsv, const char* outDir, const char* outRoot) {
  try { ep_rapidity_cf_plots::Run(binsTsv, metaTsv, panelsTsv, summaryTsv, outDir, outRoot); }
  catch (const std::exception& error) { std::cerr << "MAKE_CF_PLOTS_FAILED: " << error.what() << std::endl; gSystem->Exit(2); }
  catch (...) { std::cerr << "MAKE_CF_PLOTS_FAILED: unknown exception" << std::endl; gSystem->Exit(2); }
}

#endif // __CINT__
