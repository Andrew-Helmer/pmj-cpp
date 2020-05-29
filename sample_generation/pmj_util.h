// Copyright 2020 Andrew Helmer.
#ifndef SAMPLE_GENERATION_PMJ_UTIL_H_
#define SAMPLE_GENERATION_PMJ_UTIL_H_

#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {
  std::vector<std::pair<int, int>> GetBalancedChoicesRandom(
      const Point samples[],
      const int dim);

  std::vector<std::pair<int, int>> GetBalancedChoicesOxPlowing(
      const Point samples[],
      const int dim);

  std::vector<std::pair<int, int>> GetBalancedChoicesHilbert(
      const Point samples[],
      const int dim);
}  // namespace pmj

#endif  // SAMPLE_GENERATION_PMJ_UTIL_H_
