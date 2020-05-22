// Copyright 2020 Andrew Helmer
#ifndef PJ_H_
#define PJ_H_

#include <utility>
#include <vector>

// Progressive jittered samples shouldn't really be used, it's more just a
// learning example.
std::vector<std::pair<float, float>> prog_jittered_samples(
    const int num_samples);

#endif  // PJ_H_
