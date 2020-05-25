// Copyright 2020 Andrew Helmer
#ifndef PJ_H_
#define PJ_H_

#include <memory>
#include <utility>
#include <vector>

#include "util.h"

namespace pmj {

// Progressive jittered samples shouldn't really be used, it's more just a
// learning example.
std::unique_ptr<std::vector<Point>> GetProgJitteredSamples(
    const int num_samples);

}  // namespace pmj

#endif  // PJ_H_
