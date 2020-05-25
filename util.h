// Copyright 2020 Andrew Helmer
#ifndef UTIL_H_
#define UTIL_H_

#include <random>
#include <utility>
#include <vector>

namespace pmj {

typedef struct {
  float x;
  float y;
} Point;

// Get a random number generator. All of the functions in this file are
// thread-safe.
std::default_random_engine& get_rand_gen();

// Gets a random float between any two numbers.
float uniform_rand(float min = 0.0, float max = 1.0);

// Gets a 2D random value, based off uniform_real_distribution.
std::pair<float, float> random_sample(
    float min_x, float max_x, float min_y, float max_y);

}  // namespace pmj

#endif  // UTIL_H_
