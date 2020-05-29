// Copyright 2020 Andrew Helmer
#ifndef SAMPLE_GENERATION_PMJ02_H_
#define SAMPLE_GENERATION_PMJ02_H_

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {

// Generates progressive multi-jittered (0,2) samples WITHOUT blue noise
// properties. Takes in a number of samples.
std::unique_ptr<Point[]> GetPMJ02Samples(const int num_samples);

// Generates progressive multi-jittered (0,2) samples with blue noise
// properties.
std::unique_ptr<Point[]> GetPMJ02SamplesWithBlueNoise(
    const int num_samples);

/*
 * These functions are just for experimentation, but likely not useful for
 * real purposes, since they perform worse than the ones above.
 */

// Generates progressive multi-jittered (0,2) samples, but instead of
// ensuring (0,2) subsequences, it chooses subquadrants randomly between
// even and odd powers of two.
std::unique_ptr<Point[]> GetPMJ02SamplesNoBalance(const int num_samples);

// Generates progressive multi-jittered (0,2) samples, but instead of
// ensuring (0,2) subsequences, it chooses subquadrants using the ox-plowing
// technique in Christensen et al.
std::unique_ptr<Point[]> GetPMJ02SamplesOxPlowing(const int num_samples);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PMJ02_H_
