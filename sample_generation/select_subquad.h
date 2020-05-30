// Copyright 2020 Andrew Helmer.
#ifndef SAMPLE_GENERATION_BALANCE_UTIL_H_
#define SAMPLE_GENERATION_BALANCE_UTIL_H_

#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {
  typedef std::vector<std::pair<int, int>> (*subquad_fn)(
      const Point samples[], const int dim);
  /*
   * Pick which subquadrants to use randomly.
   */
  std::vector<std::pair<int, int>> GetSubQuadrantsRandomly(
      const Point samples[],
      const int dim);

  /*
   * This will randomly choose once to swap X or swap Y, but will then make that
   * same choice for EVERY value.
   */
  std::vector<std::pair<int, int>> GetSubQuadrantsConsistently(
    const Point samples[],
    const int dim);

  /*
   * Pick which subquadrants to use, using the ox-plowing technique from
   * Christensen et al.
   */
  std::vector<std::pair<int, int>> GetSubQuadrantsOxPlowing(
      const Point samples[],
      const int dim);
}  // namespace pmj

#endif  // SAMPLE_GENERATION_BALANCE_UTIL_H_
