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

// this increases RDF's verbosity level as long as the `verbosity` variable is
// in scope auto verbosity =
// ROOT::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(),
// ROOT::ELogLevel::kInfo);

struct Params {
  std::vector<float> func_params;
  std::vector<float> norm_params;
  std::vector<float> spline_params;
};

std::vector<Params> getRandomParams(int n) {
  TRandom3 rng;
  std::vector<Params> random_params;
  random_params.reserve(n);
  for (int i = 0; i < n; ++i) {
    Params params;
    params.func_params = {static_cast<float>(rng.Gaus(0, 0.1)),
                          static_cast<float>(rng.Gaus(0, 0.1)),
                          static_cast<float>((rng.Gaus(0, 0.2)))};
    params.norm_params = {static_cast<float>(rng.Gaus(1, 0.11)),
                          static_cast<float>(rng.Gaus(1, 0.18)),
                          static_cast<float>((rng.Gaus(1, 0.4)))};
    params.spline_params = {static_cast<float>(rng.Gaus(1, 0.3))};
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

std::vector<float> getSplineBinning(char const *filename) {
  TFile file(filename);

  auto Hist3D = file.Get<TH3F>("dev_tmp.0.0");
  std::vector<float> bins_edges;
  for (int i = 1; i <= Hist3D->GetNbinsX() + 1; ++i) {
    bins_edges.push_back(Hist3D->GetXaxis()->GetBinLowEdge(i));
  }
  return bins_edges;
}

void run_rdf_charlotte_equivalent(
    const char *dataset_name, const char *dataset_file, const Params &params,
    const std::vector<std::vector<TSpline3 *>> &splines_copies,
    const std::vector<float> &bin_edges) {

  const std::vector<std::string> cache_columns{"Enu_true", "ELep", "Q2"};

  ROOT::RDataFrame root{dataset_name, dataset_file};
  ROOT::RDF::RNode df = root.Cache<float, float, float>(cache_columns);
  df = df.Define("RecoEnu", [](float Enu_true) -> float { return Enu_true; },
                 {"Enu_true"}); // create RecoEnu columns as copy of Enu_true

  df = df.Define("ELep_shift",
                 [&params](float RecoEnu, float ELep) -> float {
                   const auto &p = params.func_params;
                   return RecoEnu + p[0] * ELep + p[1] * RecoEnu;
                 },
                 {"RecoEnu", "ELep"});

  df = df.Filter([](float Enu_true) { return Enu_true > 0 && Enu_true < 4; },
                 {"Enu_true"}, "Enu cut");

  df = df.Define("norm_weight",
                 [&params](float Q2) -> float {
                   if (Q2 < 0.25)
                     return 1.0;
                   else if (Q2 < 0.5)
                     return params.norm_params[0];
                   else if (Q2 < 2.0)
                     return params.norm_params[1];
                   else
                     return params.norm_params[2];
                 },
                 {"Q2"});

  for (auto i = 0; i < splines_copies.size(); i++) {
    const auto &splines = splines_copies[i];
    df = df.Define("spline_weight_" + std::to_string(i),
                   [&](float TrueNeutrinoEnergy) -> float {
                     auto it =
                         std::upper_bound(bin_edges.begin(), bin_edges.end(),
                                          TrueNeutrinoEnergy);
                     int bin = std::distance(bin_edges.begin(), it) - 1;

                     if (bin < 0 || bin >= (int)splines.size())
                       return 1.0; // under/overflow policy
                     return splines[bin]->Eval(params.spline_params[0]);
                   },
                   {"Enu_true"});
  }

  df = df.Define(
      "evt_weight",
      [](float norm_weight, float spline_weight_0, float spline_weight_1,
         float spline_weight_2, float spline_weight_3, float spline_weight_4,
         float spline_weight_5, float spline_weight_6, float spline_weight_7,
         float spline_weight_8, float spline_weight_9, float spline_weight_10,
         float spline_weight_11, float spline_weight_12, float spline_weight_13,
         float spline_weight_14, float spline_weight_15, float spline_weight_16,
         float spline_weight_17, float spline_weight_18, float spline_weight_19,
         float spline_weight_20, float spline_weight_21, float spline_weight_22,
         float spline_weight_23, float spline_weight_24, float spline_weight_25,
         float spline_weight_26, float spline_weight_27, float spline_weight_28,
         float spline_weight_29, float spline_weight_30, float spline_weight_31,
         float spline_weight_32, float spline_weight_33, float spline_weight_34,
         float spline_weight_35, float spline_weight_36, float spline_weight_37,
         float spline_weight_38, float spline_weight_39, float spline_weight_40,
         float spline_weight_41, float spline_weight_42, float spline_weight_43,
         float spline_weight_44, float spline_weight_45, float spline_weight_46,
         float spline_weight_47, float spline_weight_48, float spline_weight_49,
         float spline_weight_50, float spline_weight_51, float spline_weight_52,
         float spline_weight_53, float spline_weight_54, float spline_weight_55,
         float spline_weight_56, float spline_weight_57, float spline_weight_58,
         float spline_weight_59, float spline_weight_60, float spline_weight_61,
         float spline_weight_62, float spline_weight_63, float spline_weight_64,
         float spline_weight_65, float spline_weight_66, float spline_weight_67,
         float spline_weight_68, float spline_weight_69, float spline_weight_70,
         float spline_weight_71, float spline_weight_72, float spline_weight_73,
         float spline_weight_74, float spline_weight_75, float spline_weight_76,
         float spline_weight_77, float spline_weight_78, float spline_weight_79,
         float spline_weight_80, float spline_weight_81, float spline_weight_82,
         float spline_weight_83, float spline_weight_84, float spline_weight_85,
         float spline_weight_86, float spline_weight_87, float spline_weight_88,
         float spline_weight_89, float spline_weight_90, float spline_weight_91,
         float spline_weight_92, float spline_weight_93, float spline_weight_94,
         float spline_weight_95, float spline_weight_96, float spline_weight_97,
         float spline_weight_98, float spline_weight_99) {
        return norm_weight * spline_weight_0 * spline_weight_1 *
               spline_weight_2 * spline_weight_3 * spline_weight_4 *
               spline_weight_5 * spline_weight_6 * spline_weight_7 *
               spline_weight_8 * spline_weight_9 * spline_weight_10 *
               spline_weight_11 * spline_weight_12 * spline_weight_13 *
               spline_weight_14 * spline_weight_15 * spline_weight_16 *
               spline_weight_17 * spline_weight_18 * spline_weight_19 *
               spline_weight_20 * spline_weight_21 * spline_weight_22 *
               spline_weight_23 * spline_weight_24 * spline_weight_25 *
               spline_weight_26 * spline_weight_27 * spline_weight_28 *
               spline_weight_29 * spline_weight_30 * spline_weight_31 *
               spline_weight_32 * spline_weight_33 * spline_weight_34 *
               spline_weight_35 * spline_weight_36 * spline_weight_37 *
               spline_weight_38 * spline_weight_39 * spline_weight_40 *
               spline_weight_41 * spline_weight_42 * spline_weight_43 *
               spline_weight_44 * spline_weight_45 * spline_weight_46 *
               spline_weight_47 * spline_weight_48 * spline_weight_49 *
               spline_weight_50 * spline_weight_51 * spline_weight_52 *
               spline_weight_53 * spline_weight_54 * spline_weight_55 *
               spline_weight_56 * spline_weight_57 * spline_weight_58 *
               spline_weight_59 * spline_weight_60 * spline_weight_61 *
               spline_weight_62 * spline_weight_63 * spline_weight_64 *
               spline_weight_65 * spline_weight_66 * spline_weight_67 *
               spline_weight_68 * spline_weight_69 * spline_weight_70 *
               spline_weight_71 * spline_weight_72 * spline_weight_73 *
               spline_weight_74 * spline_weight_75 * spline_weight_76 *
               spline_weight_77 * spline_weight_78 * spline_weight_79 *
               spline_weight_80 * spline_weight_81 * spline_weight_82 *
               spline_weight_83 * spline_weight_84 * spline_weight_85 *
               spline_weight_86 * spline_weight_87 * spline_weight_88 *
               spline_weight_89 * spline_weight_90 * spline_weight_91 *
               spline_weight_92 * spline_weight_93 * spline_weight_94 *
               spline_weight_95 * spline_weight_96 * spline_weight_97 *
               spline_weight_98 * spline_weight_99;
      },
      {"norm_weight",      "spline_weight_0",  "spline_weight_1",
       "spline_weight_2",  "spline_weight_3",  "spline_weight_4",
       "spline_weight_5",  "spline_weight_6",  "spline_weight_7",
       "spline_weight_8",  "spline_weight_9",  "spline_weight_10",
       "spline_weight_11", "spline_weight_12", "spline_weight_13",
       "spline_weight_14", "spline_weight_15", "spline_weight_16",
       "spline_weight_17", "spline_weight_18", "spline_weight_19",
       "spline_weight_20", "spline_weight_21", "spline_weight_22",
       "spline_weight_23", "spline_weight_24", "spline_weight_25",
       "spline_weight_26", "spline_weight_27", "spline_weight_28",
       "spline_weight_29", "spline_weight_30", "spline_weight_31",
       "spline_weight_32", "spline_weight_33", "spline_weight_34",
       "spline_weight_35", "spline_weight_36", "spline_weight_37",
       "spline_weight_38", "spline_weight_39", "spline_weight_40",
       "spline_weight_41", "spline_weight_42", "spline_weight_43",
       "spline_weight_44", "spline_weight_45", "spline_weight_46",
       "spline_weight_47", "spline_weight_48", "spline_weight_49",
       "spline_weight_50", "spline_weight_51", "spline_weight_52",
       "spline_weight_53", "spline_weight_54", "spline_weight_55",
       "spline_weight_56", "spline_weight_57", "spline_weight_58",
       "spline_weight_59", "spline_weight_60", "spline_weight_61",
       "spline_weight_62", "spline_weight_63", "spline_weight_64",
       "spline_weight_65", "spline_weight_66", "spline_weight_67",
       "spline_weight_68", "spline_weight_69", "spline_weight_70",
       "spline_weight_71", "spline_weight_72", "spline_weight_73",
       "spline_weight_74", "spline_weight_75", "spline_weight_76",
       "spline_weight_77", "spline_weight_78", "spline_weight_79",
       "spline_weight_80", "spline_weight_81", "spline_weight_82",
       "spline_weight_83", "spline_weight_84", "spline_weight_85",
       "spline_weight_86", "spline_weight_87", "spline_weight_88",
       "spline_weight_89", "spline_weight_90", "spline_weight_91",
       "spline_weight_92", "spline_weight_93", "spline_weight_94",
       "spline_weight_95", "spline_weight_96", "spline_weight_97",
       "spline_weight_98", "spline_weight_99"});

  std::vector<float> bins = {0.,   0.5, 1.,   1.25, 1.5,  1.75, 2., 2.25, 2.5,
                             2.75, 3.,  3.25, 3.5,  3.75, 4.,   5., 6.,   10.};
  int nbins = bins.size() - 1;
  auto h = df.Histo1D<float, float>(
      {"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()}, "ELep_shift",
      "evt_weight");

  h->GetEntries();
}

void run_rdf(const char *dataset_name, const char *dataset_file,
             const Params &params,
             const std::vector<std::vector<TSpline3 *>> &splines_copies,
             const std::vector<float> &bin_edges) {
  // With respect to run_rdf_charlotte_equivalent, this function uses the
  // following optimisations:
  // - Reduce the number of Define nodes by computing the event weight inline in
  // one lambda instead of splitting the spline_weight computations

  const std::vector<std::string> cache_columns{"Enu_true", "ELep", "Q2"};

  ROOT::RDataFrame root{dataset_name, dataset_file};
  ROOT::RDF::RNode df = root.Cache<float, float, float>(cache_columns);
  df = df.Define("RecoEnu", [](float Enu_true) -> float { return Enu_true; },
                 {"Enu_true"}); // create RecoEnu columns as copy of Enu_true

  df = df.Define("ELep_shift",
                 [&params](float RecoEnu, float ELep) -> float {
                   const auto &p = params.func_params;
                   return RecoEnu + p[0] * ELep + p[1] * RecoEnu;
                 },
                 {"RecoEnu", "ELep"});

  df = df.Filter([](float Enu_true) { return Enu_true > 0 && Enu_true < 4; },
                 {"Enu_true"}, "Enu cut");

  df = df.Define("norm_weight",
                 [&params](float Q2) -> float {
                   if (Q2 < 0.25)
                     return 1.0;
                   else if (Q2 < 0.5)
                     return params.norm_params[0];
                   else if (Q2 < 2.0)
                     return params.norm_params[1];
                   else
                     return params.norm_params[2];
                 },
                 {"Q2"});

  df = df.Define("evt_weight",
                 [&](float norm_weight, float TrueNeutrinoEnergy) -> float {
                   auto evt_weight = norm_weight;
                   for (const auto &splines : splines_copies) {
                     auto it =
                         std::upper_bound(bin_edges.begin(), bin_edges.end(),
                                          TrueNeutrinoEnergy);
                     int bin = std::distance(bin_edges.begin(), it) - 1;

                     if (bin < 0 || bin >= (int)splines.size())
                       return 1.0; // under/overflow policy
                     evt_weight *= splines[bin]->Eval(params.spline_params[0]);
                   }
                   return evt_weight;
                 },
                 {"norm_weight", "Enu_true"});

  std::vector<float> bins = {0.,   0.5, 1.,   1.25, 1.5,  1.75, 2., 2.25, 2.5,
                             2.75, 3.,  3.25, 3.5,  3.75, 4.,   5., 6.,   10.};
  int nbins = bins.size() - 1;
  auto h = df.Histo1D<float, float>(
      {"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()}, "ELep_shift",
      "evt_weight");

  h->GetEntries();
}

void run_rntuple(const char *dataset_name, const char *dataset_file,
                 const Params &params,
                 const std::vector<std::vector<TSpline3 *>> &splines_copies,
                 const std::vector<float> &bin_edges) {

  // Create an RNTupleModel with the only three columns that will be read from
  // disk
  auto model = ROOT::RNTupleModel::Create();
  auto Enu_true = model->MakeField<float>("Enu_true");
  auto ELep = model->MakeField<float>("ELep");
  auto Q2 = model->MakeField<float>("Q2");

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

  // Open the file with the RNTuple, applying the model previously defined
  auto reader =
      ROOT::RNTupleReader::Open(std::move(model), dataset_name, dataset_file);
  std::vector<float> bins = {0.,   0.5, 1.,   1.25, 1.5,  1.75, 2., 2.25, 2.5,
                             2.75, 3.,  3.25, 3.5,  3.75, 4.,   5., 6.,   10.};
  int nbins = bins.size() - 1;
  TH1D h{"hELep", "ELep;ELep [GeV];Events", nbins, bins.data()};

  for (auto entryId : *reader) {
    reader->LoadEntry(entryId);

    // create RecoEnu as copy of Enu_true as it is done in the RDF code
    auto RecoEnu = *Enu_true;

    if (*Enu_true > 0 && *Enu_true < 4) {
      auto ELep_shift = define_ELep_shift(RecoEnu, *ELep);
      auto norm_weight = define_norm_weight(*Q2);

      float evt_weight = norm_weight;
      for (const auto &splines : splines_copies) {
        auto define_spline_weight = [&](float TrueNeutrinoEnergy) -> float {
          auto it = std::upper_bound(bin_edges.begin(), bin_edges.end(),
                                     TrueNeutrinoEnergy);
          int bin = std::distance(bin_edges.begin(), it) - 1;

          if (bin < 0 || bin >= (int)splines.size())
            return 1.0; // under/overflow policy
          return splines[bin]->Eval(params.spline_params[0]);
        };
        evt_weight *= define_spline_weight(*Enu_true);
      }

      h.Fill(ELep_shift, evt_weight);
    }
  }
}

int main() {
  // ROOT::EnableImplicitMT(8);

  auto dataset_name = "Events";
  auto dataset_file = "RNTuples/NuWro_numu_x_numu_FlatTree_Beam.root";
  auto splines_file = "BinnedSplinesTutorialInputs2D.root";

  auto splines = getSplines(splines_file);
  auto spline_binning = getSplineBinning(splines_file);

  // copy the splines so we can increase the running time/complexity
  // a real analysis will have O(100) splines
  std::vector<std::vector<TSpline3 *>> splines_copies;
  for (int i = 0; i < 100; ++i) {
    std::vector<TSpline3 *> splines_copy;
    for (auto *s : splines) {
      splines_copy.push_back(static_cast<TSpline3 *>(s->Clone()));
    }
    splines_copies.push_back(splines_copy);
  }

  // number of times to loop over the graph with different parameters,
  // equivalent to number of faked MCMC steps
  int n_trials = 1;
  auto random_params = getRandomParams(n_trials);

  auto start_rntuple = std::chrono::high_resolution_clock::now();

  run_rntuple(dataset_name, dataset_file, random_params[0], splines_copies,
              spline_binning);

  auto end_rntuple = std::chrono::high_resolution_clock::now();
  auto duration_rntuple = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_rntuple - start_rntuple);
  std::cout << "Total time (RNTuple): " << duration_rntuple.count() << " ms"
            << std::endl;

  auto start_rdf = std::chrono::high_resolution_clock::now();

  run_rdf_charlotte_equivalent(dataset_name, dataset_file, random_params[0],
                               splines_copies, spline_binning);

  auto end_rdf = std::chrono::high_resolution_clock::now();
  auto duration_rdf = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_rdf - start_rdf);
  std::cout << "Total time (RDataFrame - Charlotte equivalent): "
            << duration_rdf.count() << " ms" << std::endl;

  auto start_opt = std::chrono::high_resolution_clock::now();

  run_rdf(dataset_name, dataset_file, random_params[0], splines_copies,
          spline_binning);

  auto end_opt = std::chrono::high_resolution_clock::now();
  auto duration_opt = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_opt - start_opt);
  std::cout << "Total time (RDataFrame - Optimised): " << duration_opt.count()
            << " ms" << std::endl;
}
