// Expose only the entry point to ROOT 5 CINT. Helper templates are compiled
// by ACLiC but do not need interpreter dictionaries.
void extract_cf_inputs(const char* planPath, const char* outTsv, const char* outRoot, const char* axesTsv);

#ifndef __CINT__

#include "TArrayD.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TObject.h"
#include "TParameter.h"
#include "TSystem.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Read histogram bins only. All selections are supplied by the TSV plan;
// normalization, rebinning and correlation-function construction are separate.
namespace ep_rapidity_cf_extract {

bool rootErrorSeen = false;
std::string rootErrorMessage;

void RootErrorHandler(Int_t level, Bool_t, const char* location, const char* message) {
  if (level >= kError) {
    rootErrorSeen = true;
    std::ostringstream text;
    text << (location ? location : "ROOT") << ": " << (message ? message : "");
    rootErrorMessage = text.str();
  }
  DefaultErrorHandler(level, kFALSE, location, message);
}

struct ErrorHandlerGuard {
  ErrorHandlerFunc_t previous;
  ErrorHandlerGuard() : previous(SetErrorHandler(RootErrorHandler)) {
    rootErrorSeen = false;
    rootErrorMessage.clear();
  }
  ~ErrorHandlerGuard() { SetErrorHandler(previous); }
};

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void CheckRoot() {
  Require(!rootErrorSeen, std::string("ROOT reported an error: ") + rootErrorMessage);
}

bool Finite(double value) { return TMath::Finite(value) != 0; }

int CellCount(const TH1& histogram) {
  Long64_t cells = histogram.GetNbinsX() + 2;
  if (histogram.GetDimension() > 1) cells *= histogram.GetNbinsY() + 2;
  if (histogram.GetDimension() > 2) cells *= histogram.GetNbinsZ() + 2;
  Require(cells > 0 && cells <= std::numeric_limits<int>::max(), "Histogram cell count is outside the ROOT global-bin range");
  return static_cast<int>(cells);
}

std::string CanonicalPath(const std::string& path) {
  char resolved[PATH_MAX];
  return realpath(path.c_str(), resolved) ? std::string(resolved) : path;
}

int Integer(const std::string& text) {
  std::istringstream stream(text);
  int value = 0;
  Require(bool(stream >> value), std::string("Invalid integer: ") + text);
  stream >> std::ws;
  Require(stream.eof(), std::string("Trailing text in integer: ") + text);
  return value;
}

bool SafeDirectoryName(const std::string& name) {
  if (name.empty() || name == "." || name == "..") return false;
  for (size_t i = 0; i < name.size(); ++i) {
    const char c = name[i];
    if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
          (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.')) return false;
  }
  return true;
}

struct Request {
  std::string dataset, projection, centrality, seKey, meKey, input;
  int centMin, centMax, dimension, yFirst, yLast;
};

std::vector<Request> ReadPlan(const char* path) {
  Require(path && *path, "The plan path is required");
  std::ifstream stream(path);
  Require(bool(stream), std::string("Cannot read plan: ") + path);
  std::vector<Request> requests;
  std::set<std::string> outputs;
  std::string line;
  unsigned int lineNumber = 0;
  while (std::getline(stream, line)) {
    ++lineNumber;
    if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
    if (line.empty()) continue;
    std::vector<std::string> fields;
    size_t begin = 0;
    for (;;) {
      const size_t end = line.find('\t', begin);
      fields.push_back(line.substr(begin, end == std::string::npos ? end : end - begin));
      if (end == std::string::npos) break;
      begin = end + 1;
    }
    std::ostringstream context;
    context << "Plan line " << lineNumber;
    Require(fields.size() == 11, context.str() + " must contain 11 TSV fields");
    Request request;
    request.dataset = fields[0]; request.projection = fields[1]; request.centrality = fields[2];
    request.centMin = Integer(fields[3]); request.centMax = Integer(fields[4]);
    request.dimension = Integer(fields[5]);
    request.yFirst = Integer(fields[6]); request.yLast = Integer(fields[7]);
    request.seKey = fields[8]; request.meKey = fields[9]; request.input = CanonicalPath(fields[10]);
    Require(SafeDirectoryName(request.dataset) && SafeDirectoryName(request.projection) &&
            SafeDirectoryName(request.centrality), context.str() + " has an unsafe output-directory name");
    Require(request.centMin <= request.centMax, context.str() + " has reversed centrality bounds");
    Require(request.dimension == 2 || request.dimension == 3, context.str() + " requires dimension 2 or 3");
    Require(!request.seKey.empty() && !request.meKey.empty() && !request.input.empty(),
            context.str() + " has an empty source field");
    const std::string output = request.dataset + "/" + request.centrality + "/" + request.projection;
    Require(outputs.insert(output).second, std::string("Duplicate output selection: ") + output);
    requests.push_back(request);
  }
  Require(!stream.bad() && !requests.empty(), "The plan is empty or could not be read completely");
  return requests;
}

void ValidateAxis(const TAxis& axis, const std::string& context) {
  Require(axis.GetNbins() > 0, context + ": empty axis");
  double previous = axis.GetBinLowEdge(1);
  Require(Finite(previous), context + ": nonfinite first axis edge");
  for (int edge = 2; edge <= axis.GetNbins() + 1; ++edge) {
    const double current = axis.GetBinLowEdge(edge);
    Require(Finite(current) && current > previous, context + ": invalid axis edges");
    previous = current;
  }
}

void SameAxis(const TAxis& a, const TAxis& b, const std::string& context) {
  Require(a.GetNbins() == b.GetNbins(), context + ": SE/ME axis bin counts differ");
  for (int edge = 1; edge <= a.GetNbins() + 1; ++edge)
    Require(a.GetBinLowEdge(edge) == b.GetBinLowEdge(edge), context + ": SE/ME axis edges differ");
}

struct InputFile {
  TFile* file;
  std::string path;
  std::map<std::string, TH1*> histograms;
  std::set<std::string> recordedHistograms, recordedMetadata;
  TH1* axisReference;
  int pairBins;
  double pairMin, pairMax, he4Threshold;

  explicit InputFile(const std::string& source)
      : file(0), path(source), axisReference(0), pairBins(0), pairMin(0), pairMax(0), he4Threshold(0) {
    file = TFile::Open(path.c_str(), "READ");
    if (rootErrorSeen || !file || file->IsZombie() || file->TestBit(TFile::kRecovered)) {
      if (file) { file->Close(); delete file; file = 0; }
      CheckRoot();
      throw std::runtime_error(std::string("Invalid/recovered input ROOT: ") + path);
    }
  }

  ~InputFile() {
    for (std::map<std::string, TH1*>::iterator h = histograms.begin(); h != histograms.end(); ++h)
      delete h->second;
    if (file) { file->Close(); delete file; }
  }

  TH1* Histogram(const std::string& name, int dimension) {
    std::map<std::string, TH1*>::iterator found = histograms.find(name);
    if (found != histograms.end()) {
      Require(found->second->GetDimension() == dimension, "A cached histogram has inconsistent requested dimensions");
      return found->second;
    }
    TObject* object = file->Get(name.c_str());
    CheckRoot();
    TH1* histogram = dynamic_cast<TH1*>(object);
    Require(histogram && histogram->GetDimension() == dimension &&
            ((dimension == 2 && dynamic_cast<TH2*>(histogram)) ||
             (dimension == 3 && dynamic_cast<TH3*>(histogram))),
            std::string("Missing/wrong-dimension histogram: ") + path + ":" + name);
    histogram->SetDirectory(0);
    histograms[name] = histogram;
    ValidateAxis(*histogram->GetXaxis(), name + " X");
    ValidateAxis(*histogram->GetYaxis(), name + " Y");
    if (dimension == 3) ValidateAxis(*histogram->GetZaxis(), name + " Z");
    const int cells = CellCount(*histogram);
    Require(histogram->GetSumw2N() == 0 || histogram->GetSumw2N() == cells,
            std::string("Invalid Sumw2 array length: ") + name);
    Require(Finite(histogram->GetEntries()), std::string("Nonfinite source entries: ") + name);
    for (int cell = 0; cell < cells; ++cell) {
      const double value = histogram->GetBinContent(cell);
      const double error2 = histogram->GetSumw2N() ? histogram->GetSumw2()->At(cell) : std::fabs(value);
      Require(Finite(value) && Finite(error2) && error2 >= 0.,
              std::string("Nonfinite content or invalid variance in source: ") + name);
    }
    if (axisReference) {
      SameAxis(*axisReference->GetXaxis(), *histogram->GetXaxis(), name + " common kstar");
      const TAxis* referenceCent = axisReference->GetDimension() == 2
          ? axisReference->GetYaxis() : axisReference->GetZaxis();
      const TAxis* histogramCent = dimension == 2 ? histogram->GetYaxis() : histogram->GetZaxis();
      SameAxis(*referenceCent, *histogramCent, name + " common centrality");
    } else axisReference = histogram;
    return histogram;
  }

  void Metadata(const std::string& dataset, std::ostream& axes) {
    if (!recordedMetadata.insert(dataset).second) return;
    const char* names[] = {"he4PairRapidityLabNBins", "he4PairRapidityLabMin",
                          "he4PairRapidityLabMax", "he4RapiditySplitYLabThreshold"};
    for (int i = 0; i < 4; ++i) {
      TObject* object = file->Get(names[i]);
      CheckRoot();
      double value = 0.;
      if (i == 0) {
        TParameter<int>* parameter = dynamic_cast<TParameter<int>*>(object);
        Require(parameter != 0, std::string("Missing/wrong-type metadata: ") + names[i]);
        pairBins = parameter->GetVal(); value = pairBins;
      } else {
        TParameter<double>* parameter = dynamic_cast<TParameter<double>*>(object);
        Require(parameter != 0, std::string("Missing/wrong-type metadata: ") + names[i]);
        value = parameter->GetVal();
        if (i == 1) pairMin = value;
        if (i == 2) pairMax = value;
        if (i == 3) he4Threshold = value;
      }
      Require(Finite(value), std::string("Nonfinite metadata: ") + names[i]);
      axes << "META\t" << dataset << '\t' << path << '\t' << names[i] << '\t'
           << object->ClassName() << '\t' << value << '\n';
    }
    Require(pairBins > 0 && pairMin < pairMax, "Invalid pair-rapidity metadata");
  }

  void Record(const std::string& dataset, const std::string& key, TH1& histogram, std::ostream& axes) {
    if (!recordedHistograms.insert(dataset + "\t" + key).second) return;
    axes << "SOURCE\t" << dataset << '\t' << path << '\t' << key << '\t'
         << histogram.ClassName() << '\t' << histogram.GetDimension() << '\t'
         << histogram.GetEntries() << '\t' << histogram.GetSumw2N() << '\n';
    const TAxis* sourceAxes[] = {histogram.GetXaxis(), histogram.GetYaxis(), histogram.GetZaxis()};
    const char* labels[] = {"X", "Y", "Z"};
    for (int dimension = 0; dimension < histogram.GetDimension(); ++dimension) {
      const TAxis& axis = *sourceAxes[dimension];
      axes << "AXIS\t" << dataset << '\t' << path << '\t' << key << '\t'
           << labels[dimension] << '\t' << axis.GetNbins();
      for (int edge = 1; edge <= axis.GetNbins() + 1; ++edge) axes << '\t' << axis.GetBinLowEdge(edge);
      axes << '\n';
    }
  }

 private:
  InputFile(const InputFile&);
  InputFile& operator=(const InputFile&);
};

struct InputCache {
  std::map<std::string, InputFile*> files;
  ~InputCache() { Clear(); }
  void Clear() {
    for (std::map<std::string, InputFile*>::iterator f = files.begin(); f != files.end(); ++f) delete f->second;
    files.clear();
  }
  InputFile& Get(const std::string& path) {
    std::map<std::string, InputFile*>::iterator found = files.find(path);
    if (found != files.end()) return *found->second;
    InputFile* input = new InputFile(path);
    files[path] = input;
    return *input;
  }
};

TDirectory* Directory(TDirectory* parent, const std::string& name) {
  TDirectory* directory = parent->GetDirectory(name.c_str());
  if (!directory) directory = parent->mkdir(name.c_str());
  CheckRoot();
  Require(directory != 0, std::string("Cannot create output directory: ") + name);
  return directory;
}

void Project(TH1& source, const Request& request, const char* mode,
             TDirectory* output, std::ostream& bins, std::ostream& axes) {
  TAxis* centAxis = request.dimension == 2 ? source.GetYaxis() : source.GetZaxis();
  const int centFirst = centAxis->FindBin(request.centMin);
  const int centLast = centAxis->FindBin(request.centMax);
  Require(centFirst >= 1 && centFirst <= centLast && centLast <= centAxis->GetNbins(),
          "Centrality selection reaches a flow bin or is reversed");
  const int yFirst = request.dimension == 2 ? 1 : request.yFirst;
  const int yLast = request.dimension == 2 ? 1 : request.yLast;
  Require(request.dimension == 2 || (yFirst >= 0 && yFirst <= yLast && yLast <= source.GetNbinsY() + 1),
          "Pair-rapidity ROOT bin selection is outside [0,nY+1]");
  axes << "SELECT\t" << request.dataset << '\t' << request.projection << '\t' << request.centrality
       << '\t' << mode << '\t' << request.centMin << '\t' << request.centMax
       << '\t' << centFirst << '\t' << centLast << '\t' << request.dimension
       << '\t' << request.yFirst << '\t' << request.yLast << '\n';
  const int nX = source.GetNbinsX();
  std::vector<double> edges(nX + 1);
  for (int edge = 0; edge <= nX; ++edge) edges[edge] = source.GetXaxis()->GetBinLowEdge(edge + 1);
  const std::string name = std::string(mode) + "_raw";
  TH1D projected(name.c_str(), source.GetTitle(), nX, &edges[0]);
  projected.SetDirectory(0);
  projected.Sumw2();
  double selectedContent = 0.;
  for (int x = 0; x <= nX + 1; ++x) {
    double content = 0., variance = 0.;
    for (int cent = centFirst; cent <= centLast; ++cent) {
      for (int y = yFirst; y <= yLast; ++y) {
        const int global = request.dimension == 2 ? source.GetBin(x, cent) : source.GetBin(x, y, cent);
        const double value = source.GetBinContent(global);
        const double error2 = source.GetSumw2N() ? source.GetSumw2()->At(global) : std::fabs(value);
        Require(Finite(value) && Finite(error2) && error2 >= 0.,
                std::string("Invalid source content/variance: ") + source.GetName());
        content += value; variance += error2;
      }
    }
    Require(Finite(content) && Finite(variance) && variance >= 0., "Projection overflow/nonfinite variance");
    projected.SetBinContent(x, content);
    projected.GetSumw2()->SetAt(variance, x);
    selectedContent += content;
    const double low = source.GetXaxis()->GetBinLowEdge(x);
    const double high = source.GetXaxis()->GetBinUpEdge(x);
    Require(Finite(low) && Finite(high) && high > low, "Invalid projected flow/normal bin edges");
    bins << "BIN\t" << request.dataset << '\t' << request.projection << '\t' << request.centrality
         << '\t' << mode << '\t' << x << '\t' << low << '\t' << high
         << '\t' << content << '\t' << variance << '\n';
  }
  Require(Finite(selectedContent), "Nonfinite selected-content sum");
  // Raw source histograms contain counts. Record the selected count rather
  // than SetBinContent's artificial number of calls as the output entries.
  projected.SetEntries(selectedContent);
  output->cd();
  Require(projected.Write(name.c_str()) > 0, std::string("Cannot write projection: ") + name);
  CheckRoot();
}

struct OutputFile {
  TFile* file;
  explicit OutputFile(const char* path) : file(TFile::Open(path, "CREATE")) {
    if (rootErrorSeen || !file || file->IsZombie() || !file->IsWritable()) {
      if (file) { file->Close(); delete file; file = 0; }
      CheckRoot();
      throw std::runtime_error("Cannot create output ROOT");
    }
  }
  ~OutputFile() { if (file) { file->Close(); delete file; } }
  void Close() { if (file) { file->Flush(); file->Close(); delete file; file = 0; } CheckRoot(); }
 private:
  OutputFile(const OutputFile&);
  OutputFile& operator=(const OutputFile&);
};

void Run(const char* planPath, const char* outTsv, const char* outRoot, const char* axesTsv) {
  ErrorHandlerGuard errors;
  const Bool_t previousAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  struct AddDirectoryGuard {
    Bool_t previous;
    explicit AddDirectoryGuard(Bool_t value) : previous(value) {}
    ~AddDirectoryGuard() { TH1::AddDirectory(previous); }
  } addDirectoryGuard(previousAddDirectory);
  const std::vector<Request> requests = ReadPlan(planPath);
  const char* outputs[] = {outTsv, outRoot, axesTsv};
  std::set<std::string> outputPaths;
  for (int output = 0; output < 3; ++output) {
    Require(outputs[output] && *outputs[output], "All three output paths are required");
    const std::string path = CanonicalPath(outputs[output]);
    Require(outputPaths.insert(path).second && path != CanonicalPath(planPath), "Output paths overlap each other or the plan");
    Require(gSystem->AccessPathName(outputs[output]), std::string("Refusing to overwrite existing output: ") + outputs[output]);
    for (size_t row = 0; row < requests.size(); ++row)
      Require(path != requests[row].input, "An output path overlaps an input ROOT");
  }
  std::ofstream bins(outTsv), axes(axesTsv);
  Require(bool(bins) && bool(axes), "Cannot create output TSV files");
  bins << std::setprecision(17);
  axes << std::setprecision(17);
  OutputFile output(outRoot);
  InputCache inputs;
  for (size_t row = 0; row < requests.size(); ++row) {
    const Request& request = requests[row];
    InputFile& input = inputs.Get(request.input);
    input.Metadata(request.dataset, axes);
    TH1* se = input.Histogram(request.seKey, request.dimension);
    TH1* me = input.Histogram(request.meKey, request.dimension);
    Require(std::string(se->ClassName()) == me->ClassName(), "SE/ME histogram classes differ");
    SameAxis(*se->GetXaxis(), *me->GetXaxis(), request.projection + " X");
    SameAxis(*se->GetYaxis(), *me->GetYaxis(), request.projection + " Y");
    if (request.dimension == 3) {
      SameAxis(*se->GetZaxis(), *me->GetZaxis(), request.projection + " Z");
      Require(se->GetNbinsY() == input.pairBins && se->GetYaxis()->GetXmin() == input.pairMin &&
              se->GetYaxis()->GetXmax() == input.pairMax, "Pair-rapidity axis differs from stored metadata");
    }
    input.Record(request.dataset, request.seKey, *se, axes);
    input.Record(request.dataset, request.meKey, *me, axes);
    TDirectory* directory = Directory(Directory(Directory(output.file, request.dataset), request.centrality), request.projection);
    Project(*se, request, "SE", directory, bins, axes);
    Project(*me, request, "ME", directory, bins, axes);
    CheckRoot();
  }
  bins.flush(); axes.flush();
  Require(bool(bins) && bool(axes), "Writing TSV output failed");
  bins.close(); axes.close();
  Require(!bins.fail() && !axes.fail(), "Closing TSV output failed");
  inputs.Clear();
  CheckRoot();
  output.Close();
  std::cout << "EXTRACT_CF_INPUTS_OK requests=" << requests.size() << std::endl;
}

} // namespace ep_rapidity_cf_extract

void extract_cf_inputs(const char* planPath, const char* outTsv, const char* outRoot, const char* axesTsv) {
  try {
    ep_rapidity_cf_extract::Run(planPath, outTsv, outRoot, axesTsv);
  } catch (const std::exception& error) {
    std::cerr << "EXTRACT_CF_INPUTS_FAILED: " << error.what() << std::endl;
    gSystem->Exit(2);
  } catch (...) {
    std::cerr << "EXTRACT_CF_INPUTS_FAILED: unknown exception" << std::endl;
    gSystem->Exit(2);
  }
}

#endif // __CINT__
