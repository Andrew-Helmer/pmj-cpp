/*
 * Copyright (C) Andrew Helmer 2020.
 *
 * Licensed under MIT Open-Source License: see LICENSE. If you use this code, or
 * you generate sample sets that you use, I'd appreciate a credit in the source
 * code of your software. Just my name and/or a link to the GitHub project.
 * Thanks!
 *
 * This file implements different methods of selecting the subquadrants in
 * between odd and even powers of 4 for the PMJ and PMJ02 algorithms. Compared
 * to random, they make a big difference for the overall error!
 */
#ifndef SAMPLE_GENERATION_SELECT_SUBQUAD_H_
#define SAMPLE_GENERATION_SELECT_SUBQUAD_H_

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

#endif  // SAMPLE_GENERATION_SELECT_SUBQUAD_H_
