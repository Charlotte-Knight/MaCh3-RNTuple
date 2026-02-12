#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TSpline.h>
#include <TRandom3.h> 
#include <chrono>
#include <ROOT/RLogger.hxx>
#include <ROOT/RDFHelpers.hxx>

// this increases RDF's verbosity level as long as the `verbosity` variable is in scope
//auto verbosity = ROOT::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::ELogLevel::kDebug+10);

struct Params {
    std::vector<float> func_params;
    std::vector<float> norm_params;
    std::vector<float> spline_params;
  };

ROOT::RDF::RNode getShiftedDF(ROOT::RDF::RNode df, Params* params) {
  return df.Define("ELep_shift",
    [params](float RecoEnu, float ELep) {
      const auto& p = params->func_params;
      return RecoEnu + p[0] * ELep + p[1] * RecoEnu;
    },
    {"RecoEnu", "ELep"}
  );
}

ROOT::RDF::RNode getNormWeightedDF(ROOT::RDF::RNode df, Params* params) {
  return df.Define("norm_weight",
    [params](float Q2) -> float {
       if      (Q2 < 0.25) return 1.0;
       else if (Q2 < 0.5)  return params->norm_params[0];
       else if (Q2 < 2.0)  return params->norm_params[1];
       else                return params->norm_params[2];
    },
    {"Q2"}
  );
}

ROOT::RDF::RNode getSelectedDF(ROOT::RDF::RNode df) {
  return df.Filter("Enu_true > 0 && Enu_true < 4", "Enu cut");
}

ROOT::RDF::RResultPtr<TH1D> getHist(ROOT::RDF::RNode df) {
  std::vector<float> bins = {0., 0.5, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 5., 6., 10.};
  //std::vector<float> bins = {0., 10.};
  int nbins = bins.size() - 1;
  return df.Histo1D<float, float>(
    {"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()},
    "ELep_shift",
    "evt_weight"
  );
}

ROOT::RDF::RNode readDF(char const *filename, char const *ntuplename, std::vector<std::string> columns) {
  ROOT::RDataFrame df(ntuplename, filename);
  //auto df_cached = df;
  auto df_cached = df.Cache(columns);
  return df_cached.Define("RecoEnu", "Enu_true");
  //return df_cached.Define("RecoEnu", "Enu_true").Range(100000);
  //return df;
}

std::vector<TSpline3*> getSplines(char const *filename) {
  TFile file(filename);
  std::vector<TSpline3*> splines;
  for (int i = 0; i < 5; ++i) {
    std::string name = "dev.mysyst1.ccqe.sp." + std::to_string(i) + ".0.0";
    auto* s = file.Get<TSpline3>(name.c_str());
    splines.push_back(static_cast<TSpline3*>(s->Clone()));
  }
  return splines;
}

std::vector<float> getSplineBinning(char const *filename) {
  TFile file(filename);
  
  auto Hist3D = file.Get<TH3F>("dev_tmp.0.0");
  std::vector<float> bins_edges;
  for (int i = 1; i <= Hist3D->GetNbinsX() + 1; ++i) {
    bins_edges.push_back(Hist3D->GetXaxis()->GetBinLowEdge(i));
  }
  return bins_edges;
}

ROOT::RDF::RNode getSplineWeightedDF(ROOT::RDF::RNode df, std::vector<TSpline3*> splines, 
                                     std::vector<float> binEdges, Params* params) {
  return df.Define(
    "spline_weight",
    [splines, binEdges, params](float TrueNeutrinoEnergy) -> float {
      auto it = std::upper_bound(binEdges.begin(), binEdges.end(), TrueNeutrinoEnergy);
      int bin = std::distance(binEdges.begin(), it) - 1;

      if (bin < 0 || bin >= (int)splines.size())
        return 1.0;  // under/overflow policy

      return splines[bin]->Eval(params->spline_params[0]);
    },
    {"Enu_true"}
  );
}

std::vector<Params> getRandomParams(int n) {
  TRandom3 rng;
  std::vector<Params> random_params;
  random_params.reserve(n);
  for (int i = 0; i < n; ++i) {
    Params params; 
    params.func_params = {static_cast<float>(rng.Gaus(0, 0.1)), static_cast<float>(rng.Gaus(0, 0.1)), static_cast<float>((rng.Gaus(0, 0.2)))}; 
    params.norm_params = {static_cast<float>(rng.Gaus(1, 0.11)), static_cast<float>(rng.Gaus(1, 0.18)), static_cast<float>((rng.Gaus(1, 0.4)))}; 
    params.spline_params = {static_cast<float>(rng.Gaus(1, 0.3))};
    random_params.push_back(params);
  }
  return random_params;  
}

ROOT::RDF::RNode getReweightedDF(ROOT::RDF::RNode df, Params* params, std::vector<TSpline3*> splines, std::vector<float> binEdges) { 
  auto df_shift = getShiftedDF(df, params);
  auto df_sel = getSelectedDF(df_shift);
  auto df_norm_weighted = getNormWeightedDF(df_sel, params);
  auto df_spline_weighted = getSplineWeightedDF(df_norm_weighted, splines, binEdges, params); // add alpha as parameter if needed
  return df_spline_weighted.Define("evt_weight", "norm_weight * spline_weight");
  //return df_norm_weighted.Define("evt_weight", "norm_weight"); // ignore spline weight for now
  //return df.Define("evt_weight", []() -> float { return 1.0; }).Define("ELep_shift", "ELep"); // ignore weights for now
  //return df;
}

ROOT::RDF::RResultPtr<TH1D> getReweightedHist(ROOT::RDF::RNode df, Params* params, std::vector<TSpline3*> splines, std::vector<float> binEdges) { 
  auto df_reweighted = getReweightedDF(df, params, splines, binEdges);
  return getHist(df_reweighted); 
}

int main(int argc, char const *argv[])
{
  if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <ntuple-file-name> <spline-file-name>" << std::endl;
      return 1;
  }
  
  // ROOT::EnableImplicitMT();
  
  auto df = readDF(argv[1], "Events", {"Enu_true", "ELep", "Q2"});
  df.Describe().Print();

  auto splines = getSplines(argv[2]);
  auto spline_binning = getSplineBinning(argv[2]);

  int n_trials = 1000;
  auto random_params = getRandomParams(n_trials);

  Params* current_params = &random_params[0];

  // creates dataframe and corresponding computation graph
  auto df_reweighted = getReweightedDF(df, current_params, splines, spline_binning);
  //df_reweighted.Count().GetValue();

  // runs graph once to trigger JIT compilation 
  auto h = getHist(df_reweighted);
  h->GetEntries();

  df_reweighted.Report()->Print();

  ROOT::RDF::SaveGraph(df_reweighted, "rdf_graph.dot");

  auto start = std::chrono::high_resolution_clock::now();

  for (auto params : random_params) {
    *current_params = params;
    //df_reweighted.Count().GetValue(); // trigger computation with new parameters
    auto h = getHist(df_reweighted);
    h->GetEntries(); // force histogram to be filled
    //std::cout << "Sum: " << h->GetSum() << std::endl; // print sum to check that results change with different parameters
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Total time: " << duration.count() << " ms" << std::endl;
  std::cout << "Average time per trial: " << duration.count() / (n_trials+0.0) << " ms" << std::endl;

  return 0;
}
