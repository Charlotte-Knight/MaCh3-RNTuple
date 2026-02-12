#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TSpline.h>

struct Params {
    std::vector<double> func_params;
    std::vector<double> norm_params;
  };

ROOT::RDF::RNode getShiftedDF(ROOT::RDF::RNode df, std::vector<double> func_params) {
  return df
    .Define("ELep_shift",
      [=](float RecoEnu, float ELep) {
        return RecoEnu + func_params[0] * ELep + func_params[1] * RecoEnu;
      },
      {"RecoEnu", "ELep"}
    );
}

ROOT::RDF::RNode getNormWeightedDF(ROOT::RDF::RNode df, std::vector<double> norm_params) {
  return df.Define(
    "norm_weight",
    [=](double Q2) {
       if (Q2 < 0.25) return 1.0;
       else if (Q2 < 0.5) return norm_params[0];
       else if (Q2 < 2.0) return norm_params[1];
       else return norm_params[2];
    },
    {"Q2"}
  );
}

ROOT::RDF::RNode getSelectedDF(ROOT::RDF::RNode df) {
  return df.Filter("Enu_true > 0 && Enu_true < 4");
}

ROOT::RDF::RResultPtr<TH1D> getHist(ROOT::RDF::RNode df) {
  std::vector<double> bins = {
    0.,  0.5,  1.,  1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 5., 6., 10.
  };
  int nbins = bins.size() - 1;

  auto h = df.Histo1D(
    {"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()},
    "ELep_shift",
    "evt_weight"
  );

  return h;
}

// ROOT::RDF::RNode getReweightedDF(ROOT::RDF::RNode df, Params params, ) {
//   auto df_shift = getShiftedDF(df, params.func_params);
//   auto df_sel = getSelectedDF(df_shift);
//   auto df_weighted = getNormWeightedDF(df_sel, params.norm_params);
//   return df_weighted;
// }

ROOT::RDF::RNode readDF(char const *filename, char const *ntuplename) {
  ROOT::RDataFrame df(ntuplename, filename);
  return df.Define("RecoEnu", "Enu_true"); // add extra variables
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

std::vector<double> getSplineBinning(char const *filename) {
  TFile file(filename);
  
  auto Hist3D = file.Get<TH3F>("dev_tmp.0.0");
  std::vector<double> bins_edges;
  for (int i = 1; i <= Hist3D->GetNbinsX() + 1; ++i) {
    bins_edges.push_back(Hist3D->GetXaxis()->GetBinLowEdge(i));
  }
  return bins_edges;
}

ROOT::RDF::RNode getSplineWeightedDF(ROOT::RDF::RNode df, std::vector<TSpline3*> splines, std::vector<double> binEdges, double alpha) {
  return df.Define(
    "spline_weight",
    [splines, binEdges, alpha](double TrueNeutrinoEnergy) {
      auto it = std::upper_bound(binEdges.begin(), binEdges.end(), TrueNeutrinoEnergy);
      int bin = std::distance(binEdges.begin(), it) - 1;

      if (bin < 0 || bin >= (int)splines.size())
        return 1.0;  // under/overflow policy

      return splines[bin]->Eval(alpha);
    },
    {"Enu_true"}
  );
}

int main(int argc, char const *argv[])
{
  if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <ntuple-file-name> <spline-file-name>" << std::endl;
      return 1;
  }
  
  //ROOT::EnableImplicitMT();

  auto df = readDF(argv[1], "Events");
  df.Describe().Print();

  auto splines = getSplines(argv[2]);
  auto spline_binning = getSplineBinning(argv[2]);

  Params params;
  params.func_params = {0.1, 0.2};
  params.norm_params = {1.11, 1.18, 1.40};

  auto df_shift = getShiftedDF(df, params.func_params);
  auto df_sel = getSelectedDF(df_shift);
  auto df_norm_weighted = getNormWeightedDF(df_sel, params.norm_params);
  auto df_spline_weighted = getSplineWeightedDF(df_norm_weighted, splines, spline_binning, 1.0);

  auto df_reweighted = df_spline_weighted.Define("evt_weight", "norm_weight * spline_weight");

  df_reweighted.Display({"ELep", "RecoEnu", "ELep_shift", "Q2", "norm_weight", "spline_weight"}, 20)->Print();

  auto h = getHist(df_reweighted);

  // df_reweighted.Display({"ELep", "RecoEnu", "ELep_shift", "Q2", "evt_weight"}, 20)->Print();

  // Draw the histogram
  TCanvas c("c", "RNTuple histogram", 800, 600);
  h->Draw();
  c.SaveAs("ELep.png");
  // draw log scale
  c.SetLogy();
  c.SaveAs("ELep_log.png");

  return 0;
}
