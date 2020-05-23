// Copyright 2020 Andrew Helmer
#include "util.h"

#include <random>
#include <signal.h>
#include <utility>
#include <vector>

namespace pmj {

std::default_random_engine& get_rand_gen() {
  thread_local static std::random_device r;
  thread_local static std::default_random_engine gen(r());

  return gen;
}

float uniform_rand(float min, float max) {
  thread_local static std::uniform_real_distribution<float> uniform;

  std::uniform_real_distribution<float>::param_type param(min, max);

  float val = uniform(get_rand_gen(), param);
  if (val == max) {
    // It's insane that this is a compiler bug we need to handle,
    // see notes here:
    // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    val = std::nextafter(val, 0.0f);
  }

  return val;
}
float uniform_rand() {
  return uniform_rand(0.0, 1.0);
}

std::pair<float, float> random_sample(
    float min_x, float max_x, float min_y, float max_y) {
  return {uniform_rand(min_x, max_x), uniform_rand(min_y, max_y)};
}

}  // namespace pmj

