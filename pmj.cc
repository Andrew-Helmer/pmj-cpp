// Copyright 2020 Andrew Helmer
#include "pmj.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include "util.h"

namespace {
// Set the strata size, loop over the samples, put them in their strata.
void divide_strata(const int n,
                   const std::vector<std::pair<float, float>>& samples,
                   std::vector<bool>* x_strata,
                   std::vector<bool>* y_strata) {
  x_strata->resize(n*4);
  y_strata->resize(n*4);
  std::fill(x_strata->begin(), x_strata->end(), false);
  std::fill(y_strata->begin(), y_strata->end(), false);
  for (int i = 0; i < n; i++) {
    const auto& sample = samples[i];

    int x_strata_pos = sample.first * n * 4;
    int y_strata_pos = sample.second * n * 4;
    (*x_strata)[x_strata_pos] = true;
    (*y_strata)[y_strata_pos] = true;
  }
}

float get_1d_strata_sample(const int pos,
                           const float grid_size,
                           std::vector<bool>* strata) {
  while (true) {
    float val = uniform_rand(pos*grid_size, (pos+1)*grid_size);
    int strata_pos = val*strata->size();
    if (!(*strata)[strata_pos]) {
      (*strata)[strata_pos] = true;
      return val;
    }
  }
}

std::pair<float, float> get_sample_strata(const int x_pos,
                                          const int y_pos,
                                          const float grid_size,
                                          std::vector<bool>* x_strata,
                                          std::vector<bool>* y_strata) {
  return {get_1d_strata_sample(x_pos, grid_size, x_strata),
          get_1d_strata_sample(y_pos, grid_size, y_strata)};
}

std::pair<float, float> get_diag_sample_strata(const int x_pos,
                                               const int y_pos,
                                               const float grid_size,
                                               std::vector<bool>* x_strata,
                                               std::vector<bool>* y_strata) {
  const int diag_x_pos = x_pos ^ 1;
  const int diag_y_pos = y_pos ^ 1;

  return get_sample_strata(
      diag_x_pos, diag_y_pos, grid_size, x_strata, y_strata);
}

void get_remaining_samples(const int n,
                           const int dim,
                           const float grid_size,
                           std::vector<bool>* x_strata,
                           std::vector<bool>* y_strata,
                           std::vector<std::pair<float, float>>* samples) {
  const int num_samples = samples->size();

  // The x_balance is negative if there are more left diagonal samples than
  // right, and positive if there are more to the right. The y is similar but
  // up / down.
  int x_balance = 0;
  int y_balance = 0;
  for (int i = 0; i < n && 2*n+i < num_samples; i++) {
    const auto& sample = (*samples)[i];

    int x_pos = sample.first * dim;
    int y_pos = sample.second * dim;

    // The idea with the balance is, if we've mostly picked left cells within
    // subquadrants, we want to pick a right cell, and vice versa, and the same
    // for top and bottom. If we've netted 2 left cells and one top cell, we
    // care more about picking a right cell than a bottom cell.
    bool balance_x = abs(x_balance) > abs(y_balance);
    if (abs(x_balance) == abs(y_balance)) {
      // If they're equal, we randomly pick one to balance. This is better than
      // randomly moving in one direction, because sometimes both balances
      // equally dictate the same move.
      balance_x = uniform_rand() < 0.5;
    }
    if (balance_x) {
      bool balance_to_right = x_balance < 0;
      if (x_pos ^ balance_to_right) x_pos = x_pos ^ 1;
      else
        y_pos = y_pos ^ 1;
    } else {
      // The zeroth row is the "top" of the grid, so a higher number means
      // down.
      bool balance_to_down = y_balance < 0;
      if (y_pos ^ balance_to_down) y_pos = y_pos ^ 1;
      else
        x_pos = x_pos ^ 1;
    }

    (*samples)[2*n+i] = get_sample_strata(
        x_pos, y_pos, grid_size, x_strata, y_strata);

    // Update balances.
    if (x_pos & 1) x_balance++;
    else
      x_balance--;
    if (y_pos & 1) y_balance++;
    else
      y_balance--;

    // Get the one diagonally opposite to the one we just got.
    if (3*n+i >= num_samples) {
      continue;
    }
    (*samples)[3*n+i] = get_diag_sample_strata(
        x_pos, y_pos, grid_size, x_strata, y_strata);
  }
}

}  // namespace

std::vector<std::pair<float, float>> prog_mj_samples(
    const int num_samples) {
  std::vector<std::pair<float, float>> samples(num_samples);
  std::vector<bool> x_strata;
  std::vector<bool> y_strata;

  x_strata.reserve(num_samples);
  y_strata.reserve(num_samples);

  // Generate first sample randomly
  samples[0] = random_sample(0, 1, 0, 1);

  int n = 1;
  int dim = 2;
  float grid_size = 0.5;
  while (n < num_samples) {
    divide_strata(n, samples, &x_strata, &y_strata);

    // For every sample, we generate the diagonally opposite one at the current
    // grid level.
    for (int i = 0; i < n && n+i < num_samples; i++) {
      const auto& sample = samples[i];

      int x_pos = sample.first * dim;
      int y_pos = sample.second * dim;

      samples[n+i] = get_diag_sample_strata(
          x_pos, y_pos, grid_size, &x_strata, &y_strata);

      if (n+i >= num_samples) {
        break;
      }
    }

    get_remaining_samples(n, dim, grid_size, &x_strata, &y_strata, &samples);

    n *= 4;
    dim *= 2;
    grid_size *= 0.5;
  }

  return samples;
}
