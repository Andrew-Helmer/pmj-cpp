// Copyright 2020 Andrew Helmer
#include "sample_generation/util.h"

#include <random>
#include <signal.h>
#include <utility>
#include <vector>

namespace pmj {

namespace {
  std::default_random_engine& GetRandGen() {
    // THIS IS NOT THREAD-SAFE! These would need to be thread_local.
    static std::random_device r;
    static std::default_random_engine gen(r());

    return gen;
  }
}

double UniformRand(double min, double max) {
  // THIS IS NOT THREAD-SAFE! These would need to be thread_local.
  static std::uniform_real_distribution<double> uniform;

  std::uniform_real_distribution<double>::param_type param(min, max);

  return uniform(GetRandGen(), param);;
}

}  // namespace pmj

