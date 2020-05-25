// Copyright 2020 Andrew Helmer
#ifndef SAMPLE_GENERATION_PJ_H_
#define SAMPLE_GENERATION_PJ_H_

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {

// Progressive jittered samples shouldn't really be used, it's more just a
// learning example.
std::unique_ptr<Point[]> GetProgJitteredSamples(
    const int num_samples);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PJ_H_
