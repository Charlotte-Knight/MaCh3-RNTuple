#include <TRandom3.h> 
#include <chrono>
#include <TSystem.h>
#include <iostream>

int main() {
  int n_trials = 1000;
  int n_events = 50000;

  std::vector<std::vector<float>> ELep;
  ELep.reserve(1000);

  TRandom3 rand(0);

  for (int i = 0; i < n_trials; ++i) {
    std::vector<float> ELep_trial;
    ELep_trial.reserve(n_events);
    for (int j = 0; j < n_events; ++j) {
      ELep_trial.push_back(rand.Uniform(0., 10.));
    }
    ELep.push_back(std::move(ELep_trial));
  }

  std::vector<float> integrals;
  integrals.reserve(n_trials);

  auto start = std::chrono::high_resolution_clock::now();

  for (const auto& ELep_trial : ELep) {
    float integral = 0.;
    for (const auto& el : ELep_trial) {
      integral += el;
    }
    integrals.push_back(integral);
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
  std::cout << "Average time per trial: " << duration.count() / (n_trials+0.0) << " milliseconds" << std::endl;

  for (size_t i = 0; i < std::min(integrals.size(), size_t(100)); ++i) {
    std::cout << integrals[i] << std::endl;
  }

  return 0;
}