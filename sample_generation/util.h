// Copyright 2020 Andrew Helmer
#ifndef SAMPLE_GENERATION_UTIL_H_
#define SAMPLE_GENERATION_UTIL_H_

#include <random>
#include <utility>
#include <vector>

namespace pmj {

typedef struct {
  double x;
  double y;
} Point;

// Get a random number generator. All of the functions in this file are
// thread-safe.
std::default_random_engine& GetRandGen();

// Gets a random double between any two numbers.
double UniformRand(double min = 0.0, double max = 1.0);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_UTIL_H_
