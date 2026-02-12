#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TSpline.h>
#include <TRandom3.h> 
#include <chrono>
#include <ROOT/RLogger.hxx>

// this increases RDF's verbosity level as long as the `verbosity` variable is in scope
auto verbosity = ROOT::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::ELogLevel::kDebug+10);

int main(int argc, char const *argv[])
{
  if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <ntuple-file-name> <spline-file-name>" << std::endl;
      return 1;
  }
  auto df = ROOT::RDataFrame("Events", argv[1]);
  auto df_defined = df
    .Define("RecoEnu", 
      [](float Enu_true) -> float { return Enu_true; }, 
      {"Enu_true"})
    .Define("evt_weight",
      []() -> float { return 1.0; }
    );

  std::vector<float> bins = {
    0.,  0.5,  1.,  1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 5., 6., 10.
  };
  int nbins = bins.size() - 1;

  //auto h = df_defined.Histo1D<float>("RecoEnu");
  auto h = df_defined.Histo1D<float, float>(
    {"hRecoEnu", "RecoEnu;RecoEnu [GeV];Events", nbins, bins.data()},
    "RecoEnu",
    "evt_weight"
  );
  h->GetEntries(); // force histogram to be filled
}
