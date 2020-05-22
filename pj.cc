// Copyright 2020 Andrew Helmer
#include "pj.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include "util.h"

namespace {
void prog_jittered_samples_quadrant(
    const std::pair<float, float>& sample,
    const int n,
    const int i,
    const int x_pos,
    const int y_pos,
    const float grid_size,
    std::vector<std::pair<float, float>>* samples) {
  // Generate diagonally opposite.
  (*samples)[n+i] = get_diag_sample(x_pos, y_pos, grid_size);

  // Pick one of the two adjacent cells to generate new sample. This will go
  // much later in the sequence, so we might not do this if it's more samples
  // than requested.
  if (2*n+i > samples->size()) {
    return;
  }

  int new_x_pos = x_pos;
  int new_y_pos = y_pos;
  if (uniform_rand() < 0.5) {
    new_x_pos = x_pos ^ 1;
  } else {
    new_y_pos = y_pos ^ 1;
  }

  (*samples)[2*n+i] = get_sample(new_x_pos, new_y_pos, grid_size);

  // Do the diagonal of this one.
  if (3*n+i > samples->size()) {
    return;
  }

  (*samples)[3*n+i] = get_diag_sample(new_x_pos, new_y_pos, grid_size);
}
}  // namespace

std::vector<std::pair<float, float>> prog_jittered_samples(
    const int num_samples) {
  std::vector<std::pair<float, float>> samples(num_samples);

  // Generate first sample randomly.
  samples[0] = random_sample(0, 1, 0, 1);

  int n = 1;  // Number of samples in previous pass.
  int dim = 2;  // The number of subquadrants in one dimension.
  float grid_size = 0.5;  // The subquadrant size in one dimension, 1.0 / dim.
  while (n < num_samples) {
    for (int i = 0; i < n && n+i < num_samples; i++) {
      const auto& sample = samples[i];

      int x_pos = sample.first * dim;
      int y_pos = sample.second * dim;

      prog_jittered_samples_quadrant(
          sample, n, i, x_pos, y_pos, grid_size, &samples);
    }
    n *= 4;
    dim *= 2;
    grid_size *= 0.5;
  }

  return samples;
}
