// Copyright 2020 Andrew Helmer
#ifndef SAMPLE_GENERATION_UTIL_H_
#define SAMPLE_GENERATION_UTIL_H_

#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace pmj {

typedef struct {
  double x;
  double y;
} Point;

// Gets a random double between any two numbers. Thread-safe.
double UniformRand(double min = 0.0, double max = 1.0);

int UniformInt(int min, int max);

double GetNearestNeighborDistSq(const Point& sample,
                                const Point* sample_grid[],
                                const int dim);

typedef std::unique_ptr<pmj::Point[]> (*sample_f)(int);

sample_f GetSamplingFunction(const std::string& algorithm);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_UTIL_H_
