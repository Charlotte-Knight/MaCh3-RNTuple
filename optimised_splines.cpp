#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RLogger.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TSystem.h>
#include <chrono>

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>

#include "FastTSpline3Eval.h"

struct Params {
  std::vector<float> func_params;
  std::vector<float> norm_params;
  std::vector<float> spline_params;
};

struct RNTupleData {
  std::vector<float> Enu_true{};
  std::vector<float> ELep{};
  std::vector<float> Q2{};
};

std::vector<Params> getRandomParams(int n, int n_spline_systs) {
  TRandom3 rng;
  std::vector<Params> random_params;
  random_params.reserve(n);
  for (int i = 0; i < n; ++i) {
    Params params;
    params.func_params = {static_cast<float>(rng.Gaus(0, 0.1)),
                          static_cast<float>(rng.Gaus(0, 0.1))};
    params.norm_params = {static_cast<float>(rng.Gaus(1, 0.11)),
                          static_cast<float>(rng.Gaus(1, 0.18)),
                          static_cast<float>((rng.Gaus(1, 0.4)))};
    params.spline_params.reserve(n_spline_systs);
    for (int j = 0; j < n_spline_systs; ++j) {
      params.spline_params.push_back(static_cast<float>(rng.Gaus(1, 0.3)));
    }
    random_params.push_back(params);
  }
  return random_params;
}

std::vector<TSpline3 *> getSplines(char const *filename) {
  TFile file(filename);
  std::vector<TSpline3 *> splines;
  for (int i = 0; i < 5; ++i) {
    std::string name = "dev.mysyst1.ccqe.sp." + std::to_string(i) + ".0.0";
    auto *s = file.Get<TSpline3>(name.c_str());
    splines.push_back(static_cast<TSpline3 *>(s->Clone()));
  }
  return splines;
}

std::vector<std::vector<TSpline3 *>> getSplinesCopies(const std::vector<TSpline3 *> &splines, int n_copies) {
  std::vector<std::vector<TSpline3 *>> splines_copies;
  for (int i = 0; i < n_copies; ++i) {
    std::vector<TSpline3 *> splines_copy;
    for (auto *s : splines) {
      splines_copy.push_back(static_cast<TSpline3 *>(s->Clone()));
    }
    splines_copies.push_back(splines_copy);
  }
  return splines_copies;
}

std::vector<std::vector<FastTSpline3Eval>> getFastSplines(const std::vector<std::vector<TSpline3 *>> &splines_copies) {
  std::vector<std::vector<FastTSpline3Eval>> fast_splines_copies;
  for (const auto &splines_copy : splines_copies) {
    std::vector<FastTSpline3Eval> fast_splines_copy;
    for (auto *s : splines_copy) {
      fast_splines_copy.push_back(FastTSpline3Eval(*s));
    }
    fast_splines_copies.push_back(fast_splines_copy);
  }
  return fast_splines_copies;
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

int getSplineBin(float x, const std::vector<float> &bin_edges) {
  auto it = std::upper_bound(bin_edges.begin(), bin_edges.end(), x);
  return std::distance(bin_edges.begin(), it) - 1;
}

std::vector<int> getSplineBins(const RNTupleData &data, const std::vector<float> &bin_edges) {
  std::vector<int> bins;
  bins.reserve(data.Enu_true.size());
  for (float x : data.Enu_true) {
    bins.push_back(getSplineBin(x, bin_edges));
  }
  return bins;
}

void evaluateSplines(const std::vector<std::vector<FastTSpline3Eval>> &fast_splines, const Params &params) {
  for (int i = 0; i < fast_splines.size(); i++) {
    for (int j = 0; j < fast_splines[i].size(); j++) {
      fast_splines[i][j].Eval(params.spline_params[i]);
    }
  }
}

void printSplineValues(const std::vector<std::vector<FastTSpline3Eval>> &fast_splines){
  for (int i = 0; i < fast_splines.size(); i++) {
    for (int j = 0; j < fast_splines[i].size(); j++) {
      std::cout << "Spline " << i << ", segment " << j << ", value: " << fast_splines[i][j].GetCachedValue() << std::endl;
    }
  }
}

void run_vectors_fast(const RNTupleData &data, const Params &params,
                 const std::vector<std::vector<FastTSpline3Eval>> &fast_splines,
                 const std::vector<int> &spline_bins) {

  evaluateSplines(fast_splines, params);
  //printSplineValues(fast_splines);
  const int nSplines = fast_splines.size();

  auto define_ELep_shift = [&params](float reco_enu, float e_lep) -> float {
    const auto &p = params.func_params;
    return reco_enu + p[0] * e_lep + p[1] * reco_enu;
  };

  auto define_norm_weight = [&params](float q2) -> float {
    if (q2 < 0.25)
      return 1.0;
    else if (q2 < 0.5)
      return params.norm_params[0];
    else if (q2 < 2.0)
      return params.norm_params[1];
    else
      return params.norm_params[2];
  };

  std::vector<float> bins = {0.,   0.5, 1.,   1.25, 1.5,  1.75, 2., 2.25, 2.5,
                             2.75, 3.,  3.25, 3.5,  3.75, 4.,   5., 6.,   10.};
  int nbins = bins.size() - 1;
  TH1D h{"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()};

  for (decltype(data.Enu_true.size()) entry = 0; entry < data.Enu_true.size();
       entry++) {

    // create RecoEnu as copy of Enu_true as it is done in the RDF code
    auto RecoEnu = data.Enu_true[entry];

    if (data.Enu_true[entry] > 0 && data.Enu_true[entry] < 4) {
      auto ELep_shift = define_ELep_shift(RecoEnu, data.ELep[entry]);
      auto norm_weight = define_norm_weight(data.Q2[entry]);

      float evt_weight = norm_weight;
      //for (const auto &splines : fast_splines_copies) {
      for (int i = 0; i < nSplines; i++) {
          evt_weight *= fast_splines[i][spline_bins[entry]].GetCachedValue();
      };
      //std::cout << "ELep_shift: " << ELep_shift << ", evt_weight: " << evt_weight << std::endl;
      h.Fill(ELep_shift, evt_weight);
    }
  }
  double total = h.GetMean();
  //std::cout << total << std::endl; // Just to trigger the graph
}

ROOT::RDF::RNode get_rw_df(ROOT::RDF::RNode df, const Params* params,
             const std::vector<std::vector<FastTSpline3Eval>> &fast_splines,
             const std::vector<float> &spline_binning) {

  df = df.Define("spline_bin", [&spline_binning](float Enu_true) -> int {
    return getSplineBin(Enu_true, spline_binning);
  }, {"Enu_true"});

  df = df.Define("ELep_shift",
                 [params](float RecoEnu, float ELep) -> float {
                   const auto &p = params->func_params;
                   return RecoEnu + p[0] * ELep + p[1] * RecoEnu;
                 },
                 {"RecoEnu", "ELep"});

  df = df.Filter([](float Enu_true) { return Enu_true > 0 && Enu_true < 4; },
                 {"Enu_true"}, "Enu cut");

  df = df.Define("norm_weight",
                 [params](float Q2) -> float {
                   if (Q2 < 0.25)
                     return 1.0;
                   else if (Q2 < 0.5)
                     return params->norm_params[0];
                   else if (Q2 < 2.0)
                     return params->norm_params[1];
                   else
                     return params->norm_params[2];
                 },
                 {"Q2"});

  df = df.Define("evt_weight",
                 [fast_splines](float norm_weight, float TrueNeutrinoEnergy, int spline_bin) -> float {
                   auto evt_weight = norm_weight;

                   for (int i = 0; i < fast_splines.size(); i++) {
                    evt_weight *= fast_splines[i][spline_bin].GetCachedValue();
                   };
                   return evt_weight;
                 },
                 {"norm_weight", "Enu_true", "spline_bin"});

  return df;
}

void run_rdf_rw_fast(ROOT::RDF::RNode df_rw, 
  const std::vector<std::vector<FastTSpline3Eval>> &fast_splines, const Params* params) {

  evaluateSplines(fast_splines, *params);

  std::vector<float> bins = {0.,   0.5, 1.,   1.25, 1.5,  1.75, 2., 2.25, 2.5,
                             2.75, 3.,  3.25, 3.5,  3.75, 4.,   5., 6.,   10.};
  int nbins = bins.size() - 1;
  auto h = df_rw.Histo1D<float, float>(
      {"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()}, "ELep_shift",
      "evt_weight");

  double total = h->GetMean();
  //std::cout << total << std::endl; // Just to trigger the graph
}

RNTupleData create_rntuple_data(const char *dataset_name,
                                const char *dataset_file) {
  // Create an RNTupleModel with the only three columns that will be read from
  // disk
  auto model = ROOT::RNTupleModel::Create();
  auto Enu_true = model->MakeField<float>("Enu_true");
  auto ELep = model->MakeField<float>("ELep");
  auto Q2 = model->MakeField<float>("Q2");

  // Open the file with the RNTuple, applying the model previously defined
  auto reader =
      ROOT::RNTupleReader::Open(std::move(model), dataset_name, dataset_file);

  RNTupleData ret;

  int counter = 0;

  for (auto entryId : *reader) {
    reader->LoadEntry(entryId);

    ret.Enu_true.push_back(*Enu_true);
    ret.ELep.push_back(*ELep);
    ret.Q2.push_back(*Q2);

    counter ++;
    //if (counter > 10)     break;
  }

  return ret;
}

ROOT::RDF::RNode create_rdf(const char *dataset_name,
                            const char *dataset_file) {

  const std::vector<std::string> cache_columns{"Enu_true", "ELep", "Q2"};
  ROOT::RDataFrame root{dataset_name, dataset_file};
  ROOT::RDF::RNode df = root.Cache<float, float, float>(cache_columns);
  return df.Define("RecoEnu", [](float Enu_true) -> float { return Enu_true; },
                   {"Enu_true"}); // create RecoEnu columns as copy of Enu_true
}

void checkSplines(const std::vector<std::vector<FastTSpline3Eval>> &fast_splines_copies, const std::vector<std::vector<TSpline3 *>> &splines_copies) {
  for (size_t i = 0; i < splines_copies.size(); ++i) {
    for (size_t j = 0; j < splines_copies[i].size(); ++j) {
      std::vector<float> test_xs = {0.1f, 0.5f, 1.0f, 1.5f, 2.0f}; // Add more test points as needed
      for (const auto &x : test_xs) {
        float y_fast = fast_splines_copies[i][j].Eval(x);
        float y_slow = splines_copies[i][j]->Eval(x);
        if (std::abs(y_fast - y_slow) > 1e-5) {
          std::cerr << "Mismatch in spline evaluation at spline " << i << ", segment " << j
                  << ": fast = " << y_fast << ", slow = " << y_slow << std::endl;
        }
      }
    }
  }
}

int main() {
  //ROOT::EnableImplicitMT();

  auto dataset_name = "Events";
  auto dataset_file = "RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root";
  auto splines_file = "BinnedSplinesTutorialInputs2D.root";

  int n_spline_systs = 100;

  auto splines = getSplines(splines_file);
  auto spline_binning = getSplineBinning(splines_file);
  auto splines_copies = getSplinesCopies(splines, n_spline_systs);
  auto fast_splines = getFastSplines(splines_copies);

  checkSplines(fast_splines, splines_copies);

  // number of times to loop over the graph with different parameters,
  // equivalent to number of faked MCMC steps
  int n_trials = 100;
  auto random_params = getRandomParams(n_trials, n_spline_systs);

  // Warm up the data for both RDataFrame and standalone RNTuple+loop over
  // vectors

  auto df = create_rdf(dataset_name, dataset_file);
  df.Count().GetValue(); // Just to trigger the graph

  auto rntuple_data = create_rntuple_data(dataset_name, dataset_file);
  auto spline_bins = getSplineBins(rntuple_data, spline_binning);

  auto start_rntuple_fast = std::chrono::high_resolution_clock::now();

  std::cout << "Running vectors" << std::endl;
  for (const auto &params : random_params) {
    run_vectors_fast(rntuple_data, params, fast_splines, spline_bins);
  }

  auto end_rntuple_fast = std::chrono::high_resolution_clock::now();
  auto duration_rntuple_fast = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_rntuple_fast - start_rntuple_fast);
  std::cout << "Total time (RNTuple - Fast): " << duration_rntuple_fast.count() << " ms"
            << std::endl;
  std::cout << "Average time per trial (RNTuple - Fast): "
            << duration_rntuple_fast.count() / static_cast<double>(n_trials) << " ms"
            << std::endl;

  // -------

  Params* current_params = &random_params[0];

  auto df_rw = get_rw_df(df, current_params, fast_splines, spline_binning);
  run_rdf_rw_fast(df_rw, fast_splines, current_params);

  auto start_rw_df_fast = std::chrono::high_resolution_clock::now();

  std::cout << "Running dataframe" << std::endl;
  for (const auto &params : random_params) {
    *current_params = params;
    //run_rdf_fast(df, params, fast_splines, spline_binning);
    //printSplineValues(fast_splines);
    run_rdf_rw_fast(df_rw, fast_splines, current_params);
  }

  auto end_rw_df_fast = std::chrono::high_resolution_clock::now();
  auto duration_rw_df_fast = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_rw_df_fast - start_rw_df_fast);
  std::cout << "Total time (RDF - Fast): " << duration_rw_df_fast.count() << " ms"
            << std::endl;
  std::cout << "Average time per trial (RDF - Fast): "
            << duration_rw_df_fast.count() / static_cast<double>(n_trials) << " ms"
            << std::endl;

  return 0;
}