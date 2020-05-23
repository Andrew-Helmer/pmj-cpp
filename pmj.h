// Copyright 2020 Andrew Helmer
#ifndef PMJ_H_
#define PMJ_H_

#include <memory>
#include <utility>
#include <vector>

namespace pmj {

typedef struct {
  float x;
  float y;
} Sample;

// Generates progressive multi-jittered samples with blue noise properties.
// Takes in a number of samples.
std::unique_ptr<std::vector<Sample>> get_pmj_samples(const int num_samples);

}  // namespace pmj

#endif  // PMJ_H_
