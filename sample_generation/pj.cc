// Copyright 2020 Andrew Helmer
#include "sample_generation/pj.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {
namespace {

Point RandomSample(
    double min_x, double max_x, double min_y, double max_y) {
  return {UniformRand(min_x, max_x), UniformRand(min_y, max_y)};
}

Point GetSample(
    const int x_pos, const int y_pos, const double grid_size) {
  return RandomSample(x_pos*grid_size, (x_pos+1)*grid_size,
                      y_pos*grid_size, (y_pos+1)*grid_size);
}

void GenPJSamplesForQuadrant(
    const Point& sample,
    const int num_samples,
    const int n,
    const int i,
    const int x_pos,
    const int y_pos,
    const double grid_size,
    Point* samples) {
  // Generate diagonally opposite.
  samples[n+i] = GetSample(x_pos ^ 1, y_pos ^ 1, grid_size);

  // Pick one of the two adjacent cells to generate new sample. This will go
  // much later in the sequence, so we might not do this if it's more samples
  // than requested.
  if (2*n+i > num_samples) {
    return;
  }

  int new_x_pos = x_pos;
  int new_y_pos = y_pos;
  if (UniformRand() < 0.5) {
    new_x_pos = x_pos ^ 1;
  } else {
    new_y_pos = y_pos ^ 1;
  }

  samples[2*n+i] = GetSample(new_x_pos, new_y_pos, grid_size);

  // Do the diagonal of this one.
  if (3*n+i > num_samples) {
    return;
  }

  samples[3*n+i] = GetSample(new_x_pos ^ 1, new_y_pos ^ 1, grid_size);
}

}  // namespace

std::unique_ptr<Point[]> GetProgJitteredSamples(
    const int num_samples) {
  auto samples = std::make_unique<Point[]>(num_samples);

  // Generate first sample randomly.
  samples[0] = RandomSample(0, 1, 0, 1);

  int n = 1;  // Number of samples in previous pass.
  int dim = 2;  // The number of subquadrants in one dimension.
  double grid_size = 0.5;  // The subquadrant size in one dimension, 1.0 / dim.
  while (n < num_samples) {
    for (int i = 0; i < n && n+i < num_samples; i++) {
      const auto& sample = samples[i];

      int x_pos = sample.x * dim;
      int y_pos = sample.y * dim;

      GenPJSamplesForQuadrant(
          sample, num_samples, n, i, x_pos, y_pos, grid_size, samples.get());
    }
    n *= 4;
    dim *= 2;
    grid_size *= 0.5;
  }

  return samples;
}

}  // namespace pmj
