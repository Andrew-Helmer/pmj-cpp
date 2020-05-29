// Copyright 2020 Andrew Helmer
#ifndef SAMPLE_GENERATION_PMJ_H_
#define SAMPLE_GENERATION_PMJ_H_

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {

// Generates progressive multi-jittered samples WITHOUT blue noise properties.
// Takes in a number of samples.
std::unique_ptr<Point[]> GetProgMultiJitteredSamples(
    const int num_samples);

// Generates progressive multi-jittered samples with blue noise properties.
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoise(
    const int num_samples);

/*
 * These are for experimentation. Ox-Plowing is used by default, as it was in
 * Christensen et al. Hilbert seems better for smooth functions but worse for
 * discontinuous ones. Random is clearly the worst of the three.
 */
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesRandom(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesHilbert(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesOxPlowing(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoiseRandom(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoiseHilbert(
    const int num_samples);
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoiseOxPlowing(
    const int num_samples);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PMJ_H_
