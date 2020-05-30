/*
 * Copyright (C) Andrew Helmer 2020.
 *
 * Licensed under MIT Open-Source License: see LICENSE. If you use this code, or
 * you generate sample sets that you use, I'd appreciate a credit in the source
 * code of your software. Just my name and/or a link to the GitHub project.
 * Thanks!
 *
 * Generate Progressive Multi-Jittered Sequences from
 * "Progressive Multi-Jittered Sample Sequences", Christensen et al. 2018. The
 * non-best-candidate sequences generate about 2 million samples/sec for me,
 * it's not very optimized compared to Christensen's paper.
 * The best candidate sequences only do about ~230k samples/sec.
 *
 * If you're reading this code for the first time and want to understand the
 * algorithm, start with the function "GenerateSamples".
 *
 */
#ifndef SAMPLE_GENERATION_PMJ_H_
#define SAMPLE_GENERATION_PMJ_H_

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {

// Generates progressive multi-jittered samples without blue noise properties.
// Takes in a number of samples.
std::unique_ptr<Point[]> GetProgMultiJitteredSamples(
    const int num_samples);

// Generates progressive multi-jittered samples with blue noise properties, i.e.
// using best-candidate points.
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoise(
    const int num_samples);

/*
 * --------------------------------------------------------------------------
 * These are for experimentation. Ox-Plowing is used by default, as it was in
 * Christensen et al. we just have an explicit function here for clarity.
 * Random performs much worse than Ox-Plowing.
 * --------------------------------------------------------------------------
 */
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesRandom(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesOxPlowing(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoiseRandom(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoiseOxPlowing(
    const int num_samples);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PMJ_H_
