// Copyright 2020 Andrew Helmer
#include "util.h"

#include <random>
#include <signal.h>
#include <utility>
#include <vector>

namespace pmj {

std::default_random_engine& GetRandGen() {
  thread_local static std::random_device r;
  thread_local static std::default_random_engine gen(r());

  return gen;
}

double UniformRand(double min, double max) {
  thread_local static std::uniform_real_distribution<double> uniform;

  std::uniform_real_distribution<double>::param_type param(min, max);

  return uniform(GetRandGen(), param);;
}

}  // namespace pmj

