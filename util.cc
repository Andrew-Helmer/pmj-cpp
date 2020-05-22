// Copyright 2020 Andrew Helmer
#include "util.h"

#include <random>
#include <utility>
#include <vector>

std::default_random_engine& get_rand_gen() {
  thread_local static std::random_device r;
  thread_local static std::default_random_engine gen(r());

  return gen;
}

float uniform_rand(float min, float max) {
  thread_local static std::uniform_real_distribution<float> uniform;

  std::uniform_real_distribution<float>::param_type param(min, max);

  return uniform(get_rand_gen(), param);
}
float uniform_rand() {
  return uniform_rand(0.0, 1.0);
}

std::pair<float, float> random_sample(
    float min_x, float max_x, float min_y, float max_y) {
  return {uniform_rand(min_x, max_x), uniform_rand(min_y, max_y)};
}
