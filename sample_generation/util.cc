// Copyright 2020 Andrew Helmer
#include "sample_generation/util.h"

#include <random>
#include <signal.h>
#include <utility>
#include <vector>

namespace pmj {

thread_local static std::random_device r;
thread_local static std::default_random_engine gen(r());

double UniformRand(double min, double max) {
  thread_local static std::uniform_real_distribution<double> uniform;

  std::uniform_real_distribution<double>::param_type param(min, max);

  return uniform(gen, param);
}

int UniformInt(int min, int max) {
  thread_local static std::uniform_int_distribution<int> uniform;

  std::uniform_int_distribution<int>::param_type param(min, max);

  return uniform(gen, param);
}

}  // namespace pmj

