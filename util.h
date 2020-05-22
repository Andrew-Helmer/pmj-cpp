// Copyright 2020 Andrew Helmer
#ifndef UTIL_H_
#define UTIL_H_

#include <random>
#include <utility>
#include <vector>

namespace pmj {

// Get a random number generator. All of the functions in this file are
// thread-safe.
std::default_random_engine& get_rand_gen();

// Gets a random float from 0.0 to 1.0. Based on uniform_real_distribution.
float uniform_rand();

// Gets a random float between any two numbers.
float uniform_rand(float min, float max);

// Gets a 2D random value, based off uniform_real_distribution.
std::pair<float, float> random_sample(
    float min_x, float max_x, float min_y, float max_y);

}  // namespace pmj

#endif  // UTIL_H_
