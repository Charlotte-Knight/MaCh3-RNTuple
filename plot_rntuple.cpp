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
  };

ROOT::RDF::RNode getShiftedDF(ROOT::RDF::RNode df, std::vector<float> func_params) {
  return df
    .Define("ELep_shift",
      [=](float RecoEnu, float ELep) -> float {
        return RecoEnu + func_params[0] * ELep + func_params[1] * RecoEnu;
      },
      {"RecoEnu", "ELep"}
    );
}

ROOT::RDF::RNode getNormWeightedDF(ROOT::RDF::RNode df, std::vector<float> norm_params) {
  return df.Define(
    "norm_weight",
    [=](float Q2) -> float {
       if (Q2 < 0.25) return 1.0;
       else if (Q2 < 0.5) return norm_params[0];
       else if (Q2 < 2.0) return norm_params[1];
       else return norm_params[2];
    },
    {"Q2"}
  );
}

ROOT::RDF::RNode getSelectedDF(ROOT::RDF::RNode df) {
  return df.Filter(
    [](float Enu_true) -> bool { return Enu_true > 0 && Enu_true < 4; },
    {"Enu_true"}
  );
}

ROOT::RDF::RResultPtr<TH1D> getHist(ROOT::RDF::RNode df) {
  std::vector<float> bins = {
    0.,  0.5,  1.,  1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 5., 6., 10.
  };
  int nbins = bins.size() - 1;

  auto h = df.Histo1D<float, float>(
    {"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()},
    "ELep_shift",
    "evt_weight"
  );

  return h;
}

ROOT::RDF::RNode readDF(char const *filename, char const *ntuplename, std::vector<std::string> columns) {
  ROOT::RDataFrame df(ntuplename, filename);
  auto df_cached = df.Cache(columns);
  return df_cached.Define(
    "RecoEnu",
    [](float Enu_true) -> float { return Enu_true; }, 
    {"Enu_true"}); // add extra variables
}

std::vector<TSpline3*> getSplines(char const *filename) {
  TFile file(filename);
  std::vector<TSpline3*> splines;
  for (int i = 0; i < 5; ++i) {
    std::string name =
      "dev.mysyst1.ccqe.sp." + std::to_string(i) + ".0.0";
    
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

ROOT::RDF::RNode getSplineWeightedDF(ROOT::RDF::RNode df, std::vector<TSpline3*> splines, std::vector<float> binEdges, float alpha) {
  return df.Define(
    "spline_weight",
    [splines, binEdges, alpha](float TrueNeutrinoEnergy) -> float {
      auto it = std::upper_bound(binEdges.begin(), binEdges.end(), TrueNeutrinoEnergy);
      int bin = std::distance(binEdges.begin(), it) - 1;

      if (bin < 0 || bin >= (int)splines.size())
        return 1.0;  // under/overflow policy

      return splines[bin]->Eval(alpha);
    },
    {"Enu_true"}
  );
}

Params getRandomParams() {
  TRandom3 rng;
  Params params; 
  params.func_params = {static_cast<float>(rng.Gaus(0, 0.1)), static_cast<float>(rng.Gaus(0, 0.1)), static_cast<float>((rng.Gaus(0, 0.2)))}; 
  params.norm_params = {static_cast<float>(rng.Gaus(1, 0.11)), static_cast<float>(rng.Gaus(1, 0.18)), static_cast<float>((rng.Gaus(1, 0.4)))}; 
  return params;  
}

ROOT::RDF::RNode getReweightedDF(ROOT::RDF::RNode df, Params params, std::vector<TSpline3*> splines, std::vector<float> binEdges) { 
  auto df_shift = getShiftedDF(df, params.func_params);
  auto df_sel = getSelectedDF(df_shift);
  auto df_norm_weighted = getNormWeightedDF(df_sel, params.norm_params);
  auto df_spline_weighted = getSplineWeightedDF(df_norm_weighted, splines, binEdges, 1.0); 
  return df_spline_weighted.Define(
    "evt_weight", 
    [](float norm_weight, float spline_weight) -> float { return norm_weight * spline_weight; }
    ,{"norm_weight", "spline_weight"} 
  );
}

ROOT::RDF::RResultPtr<TH1D> getReweightedHist(ROOT::RDF::RNode df, Params params, std::vector<TSpline3*> splines, std::vector<float> binEdges) { 
  auto df_reweighted = getReweightedDF(df, params, splines, binEdges);
  return getHist(df_reweighted); 
}

int main(int argc, char const *argv[])
{
  if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <ntuple-file-name> <spline-file-name>" << std::endl;
      return 1;
  }
  
  ROOT::EnableImplicitMT();
  
  auto df = readDF(argv[1], "Events", {"Enu_true", "ELep", "Q2"});
  df.Describe().Print();

  auto splines = getSplines(argv[2]);
  auto spline_binning = getSplineBinning(argv[2]);

  int n_trials = 100;
  std::vector<Params> random_params;
  random_params.reserve(n_trials);
  for (int i = 0; i < n_trials; ++i) { 
    random_params.push_back(getRandomParams()); 
  }

  // force loading of data?
  auto h = getReweightedHist(df, random_params[0], splines, spline_binning);
  h->GetEntries();

  ROOT::RDF::SaveGraph(h, "rdf_graph.dot");

  auto start = std::chrono::high_resolution_clock::now();

  for (const auto& params : random_params) {
    auto h = getReweightedHist(df, params, splines, spline_binning);
    h->GetEntries(); // force histogram to be filled
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Total time: " << duration.count() << " ms" << std::endl;
  std::cout << "Average time per trial: " << duration.count() / n_trials << " ms" << std::endl;

  // df_reweighted.Display({"ELep", "RecoEnu", "ELep_shift", "Q2", "evt_weight"}, 20)->Print();

  // Draw the histogram
  // TCanvas c("c", "RNTuple histogram", 800, 600);
  // h->Draw();
  // c.SaveAs("ELep.png");
  // // draw log scale
  // c.SetLogy();
  // c.SaveAs("ELep_log.png");

  return 0;
}
