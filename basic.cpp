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

int main(int argc, char const *argv[])
{
  if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <ntuple-file-name> <spline-file-name>" << std::endl;
      return 1;
  }
  
  ROOT::RDataFrame df("Events", argv[1]);
  //auto df_cached = df.Cache<float, float, float, float>({"ELep", "Enu_true", "Q2", "CosLep"});
  //auto df_cached = df.Cache<float>({"ELep"});
  auto df_cached = df.Cache({"ELep"});
  //auto df_cached = df;

  auto h = df_cached.Histo1D<float>(
      {"hELep", "ELep;ELep [GeV];Events", 1, 0., 10.},
      "ELep" );
  h->Integral();

  int n_trials = 5000;
  std::vector<float> integrals;

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < n_trials; ++i) {
    h = df_cached.Histo1D<float>(
      {"hELep", "ELep;ELep [GeV];Events", 1, 0., 10.},
      "ELep" );
    integrals.push_back(h->Integral());
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Total time: " << duration.count() << " ms" << std::endl;
  std::cout << "Average time per trial: " << duration.count() / (n_trials+0.0) << " ms" << std::endl;

  return 0;
}