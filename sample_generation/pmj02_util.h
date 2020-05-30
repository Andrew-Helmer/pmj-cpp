/*
 * Copyright (C) Andrew Helmer 2020.
 *
 * Licensed under MIT Open-Source License: see LICENSE. If you use this code, or
 * you generate sample sets that you use, I'd appreciate a credit in the source
 * code of your software. Just my name and/or a link to the GitHub project.
 * Thanks!
 * Functionality used for PMJ(0,2) sampling.
 */
#ifndef SAMPLE_GENERATION_PMJ02_UTIL_H_
#define SAMPLE_GENERATION_PMJ02_UTIL_H_

#include <utility>
#include <vector>

namespace pmj {

// Implementation of "Efficient Generation of Points that Satisfy
// Two-Dimensional Elementary Intervals" by Matt Pharr (2019). Given a set of
// strata corresponding to all elementary (0,2) intervals, and a set of
// positions on a square grid, returns the narrowest X strata and Y strata that
// are unoccupied.
std::pair<std::vector<int>, std::vector<int>>
    GetValidStrata(const int x_pos,
                   const int y_pos,
                   const std::vector<std::vector<bool>>& strata);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PMJ02_UTIL_H_
