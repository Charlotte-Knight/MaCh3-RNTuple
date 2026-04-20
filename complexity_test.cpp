#include <ROOT/RDataFrame.hxx>
#include <chrono>

float sum(std::vector<float> vec) {
  float total = 0.0f;
  for (float v : vec) {
    total += v;
  }
  return total;
}

void run_vectors(const std::vector<float> &ELep, int n_sin) {
  std::vector<float> bins = {0.,   0.5, 1.,   1.25, 1.5,  1.75, 2., 2.25, 2.5,
                             2.75, 3.,  3.25, 3.5,  3.75, 4.,   5., 6.,   10.};
  int nbins = bins.size() - 1;
  TH1D h{"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()};

  for (float e_lep : ELep) {
    float w = 1.0;
    for (int i = 0; i < n_sin; i++) {
        w *= TMath::Sin(e_lep + i) + 1.0; // some arbitrary weight function that is more expensive to compute
    }
    h.Fill(e_lep, w);
  }
}

ROOT::RDF::RNode rw_rdf(ROOT::RDF::RNode df, int n_sin) {
  return df.Define("evt_weight", [n_sin](float ELep) -> float {
    float w = 1.0;
    for (int i = 0; i < n_sin; i++) {
        w *= TMath::Sin(ELep + i) + 1.0; // some arbitrary weight function that is more expensive to compute
    }
    return w; 
    }, {"ELep"});
}

ROOT::RDF::RResultPtr<TH1D> getHist(ROOT::RDF::RNode df) {
  std::vector<float> bins = {0.,   0.5, 1.,   1.25, 1.5,  1.75, 2., 2.25, 2.5,
                             2.75, 3.,  3.25, 3.5,  3.75, 4.,   5., 6.,   10.};
  int nbins = bins.size() - 1;
  return df.Histo1D<float, float>(
    {"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()},
    "ELep",
    "evt_weight"
  );
}

int main(int argc, char const *argv[])
{ 
  //ROOT::EnableImplicitMT(8);

  ROOT::RDataFrame df("Events", argv[1]);
  auto df_cached = df.Cache<float>({"ELep"});
  auto ELep = df_cached.Take<float>("ELep").GetValue();

  std::cout << "Number of events: " << ELep.size() << std::endl;

  int n_trials = 10;

  // -------

  std::vector<int> n_sin_values = {1, 5, 10, 50, 100};

  for (int n_sin : n_sin_values) {
    std::cout << "Running with n_sin = " << n_sin << std::endl;

    std::cout << "  Running vectors" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n_trials; ++i) {
      run_vectors(ELep, n_sin);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration_milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    auto duration_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    //std::cout << "Total time: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "    Average time per trial: " << duration_milli.count() / n_trials << " milliseconds" << std::endl;
    std::cout << "    Time per trial per event: " << duration_nano.count() / (n_trials * ELep.size()) << " nanoseconds" << std::endl;

    // -------

    auto df_reweighted = rw_rdf(df_cached, n_sin);
    auto h = getHist(df_reweighted);
    h->GetEntries(); // trigger JIT compilation

    std::cout << "  Running RDataFrame" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n_trials; ++i) {
      h = getHist(df_reweighted);
      h->GetEntries(); // trigger execution of graph
    }

    end = std::chrono::high_resolution_clock::now();
    duration_milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    duration_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    //std::cout << "Total time: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "    Average time per trial: " << duration_milli.count() / n_trials << " milliseconds" << std::endl;
    std::cout << "    Time per trial per event: " << duration_nano.count() / (n_trials * ELep.size()) << " nanoseconds" << std::endl;
  }

  return 0;
}
