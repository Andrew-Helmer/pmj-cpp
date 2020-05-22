// Copyright 2020 Andrew Helmer
#ifndef PMJ_H_
#define PMJ_H_

#include <utility>
#include <vector>

// Generates progressive multi-jittered samples with blue noise properties.
// Takes in a number of samples.
std::vector<std::pair<float, float>> prog_mj_samples(
    const int num_samples);

#endif  // PMJ_H_
