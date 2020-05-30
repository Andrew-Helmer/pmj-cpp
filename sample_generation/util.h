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
// Generates a random int in the given range. Thread-safe.
int UniformInt(int min, int max);

// Given a set of samples, a grid that points to existing samples, and the
// number of cells in one dimension of that grid, returns the candidate which
// is the furthest from all existing points.
Point GetBestCandidateOfSamples(const std::vector<Point>& candidates,
                                const Point* sample_grid[],
                                const int dim);

// Given a string like "pmj" or "pmj02bn", returns the function to generate
// samples from that algorithm.
typedef std::unique_ptr<pmj::Point[]> (*sample_fn)(int);
sample_fn GetSamplingFunction(const std::string& algorithm);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_UTIL_H_
