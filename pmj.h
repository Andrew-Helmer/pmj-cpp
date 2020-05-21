// Copyright 2020 Andrew Helmer
#ifndef PMJ_H_
#define PMJ_H_

#include <utility>
#include <vector>

// Progressive jittered samples shouldn't really be used, it's more just a
// learning example.
std::vector<std::pair<float, float>> prog_jittered_samples(
    int num_samples);
// Generates progressive multi-jittered samples. Takes in a number of samples.
// Note that it's generally recommended to use powers of 2 for the number of
// samples.
std::vector<std::pair<float, float>> prog_mj_samples(
    int num_samples);
std::vector<std::pair<float, float>> prog_mj_samples_blue_noise(
    int num_samples);

#endif  // PMJ_H_
