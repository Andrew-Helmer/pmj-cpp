// Copyright 2020 Andrew Helmer
#ifndef PMJ_H_
#define PMJ_H_

#include <memory>
#include <utility>
#include <vector>

#include "util.h"

namespace pmj {

// Generates progressive multi-jittered samples WITHOUT blue noise properties.
// Takes in a number of samples.
std::unique_ptr<std::vector<Point>> get_pmj_samples(const int num_samples);

// Generates progressive multi-jittered samples with blue noise properties.
std::unique_ptr<std::vector<Point>> get_best_candidate_pmj_samples(
    const int num_samples);

}  // namespace pmj

#endif  // PMJ_H_
