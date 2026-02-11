#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>

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
    "evt_weight",
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

ROOT::RDF::RNode getReweightedDF(ROOT::RDF::RNode df, Params params) {
  auto df_shift = getShiftedDF(df, params.func_params);
  auto df_sel = getSelectedDF(df_shift);
  auto df_weighted = getNormWeightedDF(df_sel, params.norm_params);
  return df_weighted;
}

ROOT::RDF::RNode readDF(char const *filename, char const *ntuplename) {
  ROOT::RDataFrame df(ntuplename, filename);
  return df.Define("RecoEnu", "Enu_true"); // add extra variables
}

int main(int argc, char const *argv[])
{
  if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <ntuple-file-name>" << std::endl;
      return 1;
  }
  
  ROOT::EnableImplicitMT();

  auto df = readDF(argv[1], "Events");
  df.Describe().Print();

  Params params;
  params.func_params = {0.1, 0.2};
  params.norm_params = {1.11, 1.18, 1.40};

  auto df_reweighted = getReweightedDF(df, params);
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
