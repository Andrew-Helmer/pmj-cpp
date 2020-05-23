// Copyright 2020 Andrew Helmer
#include "pmj.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

#include "util.h"

namespace pmj {
namespace {

class SampleSet {
 public:
  explicit SampleSet(const int num_samples) : num_samples(num_samples) {
    samples_ = std::make_unique<std::vector<Sample>>();
    samples_->reserve(num_samples);
    samples_->push_back({0.0, 0.0});

    x_strata.reserve(num_samples);
    x_strata = {false};

    y_strata.reserve(num_samples);
    y_strata = {false};

    sample_grid.reserve(num_samples);
    sample_grid = {nullptr};

    n_ = 1;
    dim_ = 1;
    grid_size_ = 1.0;
  }

  void create_new_sample(const int i, const int x_pos, const int y_pos);

  void subdivide_strata();

  const Sample& sample(int i) {
    return (*samples_)[i];
  }
  std::unique_ptr<std::vector<Sample>> release_samples() {
    auto samples = std::make_unique<std::vector<Sample>>();
    samples.swap(samples_);
    return samples;
  }

  const int n() { return n_; }
  const int dim() { return dim_; }
  const float grid_size() { return grid_size_; }

  std::vector<bool> x_strata;
  std::vector<bool> y_strata;
  std::vector<const Sample*> sample_grid;

  const int num_samples;

 private:
  void set_sample(const int i, const Sample& sample);

  int n_;  // Number of samples in the next pass.
  int dim_;  // Number of cells in one dimension in next pass, i.e. sqrt(n).
  float grid_size_;  // 1.0 / dim_
  std::unique_ptr<std::vector<Sample>> samples_;
};

void SampleSet::subdivide_strata() {
  const float old_n = n_;

  n_ *= 4;
  dim_ *= 2;
  grid_size_ *= 0.5;

  samples_->resize(n_);

  x_strata.resize(n_);
  y_strata.resize(n_);
  sample_grid.resize(n_);

  std::fill(sample_grid.begin(), sample_grid.end(), nullptr);
  std::fill(x_strata.begin(), x_strata.end(), false);
  std::fill(y_strata.begin(), y_strata.end(), false);
  for (int i = 0; i < old_n; i++) {
    const auto& sample = (*samples_)[i];

    x_strata[sample.x * n_] = true;
    y_strata[sample.y * n_] = true;

    int x_pos = sample.x * dim_;
    int y_pos = sample.y * dim_;
    sample_grid[y_pos*dim_ + x_pos] = &sample;
  }
}

// This generates a sample within the grid position, verifying that it doesn't
// overlap strata with any other sample.
float get_1d_strata_sample(const int pos,
                           const float grid_size,
                           const std::vector<bool>& strata) {
  while (true) {
    float val = uniform_rand(pos*grid_size, (pos+1)*grid_size);
    int strata_pos = val*strata.size();
    if (!strata[strata_pos]) {
      return val;
    }
  }
}

void SampleSet::create_new_sample(const int i,
                                  const int x_pos,
                                  const int y_pos) {
  (*samples_)[i] = {get_1d_strata_sample(x_pos, grid_size_, x_strata),
                       get_1d_strata_sample(y_pos, grid_size_, y_strata)};
}

void SampleSet::set_sample(const int i,
                           const Sample& sample) {
  (*samples_)[i] = sample;

  x_strata[sample.x * n_] = true;
  y_strata[sample.y * n_] = true;

  int x_pos = sample.x * dim_;
  int y_pos = sample.y * dim_;
  sample_grid[y_pos*dim_ + x_pos] = &(*samples_)[i];
}

void get_remaining_samples(const int n,
                           const int dim,
                           SampleSet* sample_set) {
  const int num_samples = sample_set->num_samples;

  // The x_balance is negative if there are more left diagonal samples than
  // right, and positive if there are more to the right. The y is similar but
  // up / down.
  int x_balance = 0;
  int y_balance = 0;
  for (int i = 0; i < n && 2*n+i < num_samples; i++) {
    const auto& sample = sample_set->sample(i);

    int x_pos = sample.x * dim;
    int y_pos = sample.y * dim;

    // The idea with the balance is, if we've mostly picked left cells within
    // subquadrants, we want to pick a right cell, and vice versa, and the same
    // for top and bottom. If we've netted 2 left cells and one top cell, we
    // care more about picking a right cell than a bottom cell.
    bool use_x_balance = abs(x_balance) > abs(y_balance);
    if (abs(x_balance) == abs(y_balance)) {
      // If they're equal, we randomly pick one to balance. This is better than
      // randomly moving in one direction, because sometimes both balances
      // equally dictate the same move.
      use_x_balance = uniform_rand() < 0.5;
    }
    if (use_x_balance) {
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

    sample_set->create_new_sample(2*n+i, x_pos, y_pos);

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

    sample_set->create_new_sample(3*n+i, x_pos ^ 1, y_pos ^ 1);
  }
}

}  // namespace

std::unique_ptr<std::vector<Sample>> get_pmj_samples(
    const int num_samples) {
  SampleSet sample_set(num_samples);

  // Generate first sample.
  sample_set.create_new_sample(0, 0, 0);

  int n = 1;
  while (n < num_samples) {
    sample_set.subdivide_strata();

    // For every sample, we generate the diagonally opposite one at the current
    // grid level.
    for (int i = 0; i < n && n+i < num_samples; i++) {
      const auto& sample = sample_set.sample(i);

      int x_pos = sample.x * sample_set.dim();
      int y_pos = sample.y * sample_set.dim();

      sample_set.create_new_sample(n+i, x_pos ^ 1, y_pos ^ 1);
      if (n+i >= num_samples) {
        break;
      }
    }

    get_remaining_samples(n, sample_set.dim(), &sample_set);
    n *= 4;
  }

  return sample_set.release_samples();
}

}  // namespace pmj
