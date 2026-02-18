#include <ROOT/RDataFrame.hxx>
#include <chrono>

float sum(std::vector<float> vec) {
  float total = 0.0f;
  for (float v : vec) {
    total += v;
  }
  return total;
}

int main(int argc, char const *argv[])
{  
  ROOT::RDataFrame df("Events", argv[1]);
  auto df_cached = df.Cache<float>({"ELep"});
  auto ELep = df_cached.Take<float>("ELep").GetValue();

  auto sum = df_cached.Sum<float>("ELep");

  int n_trials = 10000;
  std::vector<float> integrals;
  integrals.reserve(n_trials);

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < n_trials; ++i) {
    sum = df_cached.Sum<float>("ELep");
    integrals.push_back(sum.GetValue());
    //integrals.push_back(sum(ELep));
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "Total time: " << duration.count() << " microseconds" << std::endl;
  std::cout << "Average time per trial: " << duration.count() / n_trials << " microseconds" << std::endl;

  return 0;
}