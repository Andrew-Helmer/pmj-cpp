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

    x_strata_.reserve(num_samples);
    x_strata_ = {false};

    y_strata_.reserve(num_samples);
    y_strata_ = {false};

    sample_grid_.reserve(num_samples);
    sample_grid_ = {nullptr};

    n_ = 1;
    dim_ = 1;
    grid_size_ = 1.0;
  }

  // This generates a new sample at the current index, given the X position
  // and Y position of the subquadrant. It won't generate a new sample in an
  // existing strata.
  void create_new_sample(const int sample_index,
                         const int x_pos,
                         const int y_pos);

  void subdivide_strata();

  const Sample& sample(const int sample_index) {
    return (*samples_)[sample_index];
  }
  std::unique_ptr<std::vector<Sample>> release_samples() {
    auto samples = std::make_unique<std::vector<Sample>>();
    samples.swap(samples_);
    return samples;
  }

  const int n() { return n_; }
  const int dim() { return dim_; }
  const float grid_size() { return grid_size_; }

  const int num_samples;

 private:
  void set_sample(const int i, const Sample& sample);

  float get_nearest_neighbor_dist_sq(const float x, const float y);

  std::unique_ptr<std::vector<Sample>> samples_;

  std::vector<bool> x_strata_;
  std::vector<bool> y_strata_;

  // The sample grid is used for nearest neighbor lookups.
  std::vector<const Sample*> sample_grid_;

  int n_;  // Number of samples in the next pass.
  int dim_;  // Number of cells in one dimension in next pass, i.e. sqrt(n).
  float grid_size_;  // 1.0 / dim_

  const int n_best_candidates_ = 20;
};

void SampleSet::subdivide_strata() {
  const float old_n = n_;

  n_ *= 4;
  dim_ *= 2;
  grid_size_ *= 0.5;

  samples_->resize(n_);

  x_strata_.resize(n_);
  y_strata_.resize(n_);
  sample_grid_.resize(n_);

  std::fill(sample_grid_.begin(), sample_grid_.end(), nullptr);
  std::fill(x_strata_.begin(), x_strata_.end(), false);
  std::fill(y_strata_.begin(), y_strata_.end(), false);
  for (int i = 0; i < old_n; i++) {
    const auto& sample = (*samples_)[i];

    x_strata_[sample.x * n_] = true;
    y_strata_[sample.y * n_] = true;

    int x_pos = sample.x * dim_;
    int y_pos = sample.y * dim_;
    sample_grid_[y_pos*dim_ + x_pos] = &sample;
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

void SampleSet::create_new_sample(const int sample_index,
                                  const int x_pos,
                                  const int y_pos) {
  Sample max_sample;
  float max_dist_sq = 0.0;
  for (int i = 0; i < n_best_candidates_; i++) {
    Sample cand_sample = {get_1d_strata_sample(x_pos, grid_size_, x_strata_),
                          get_1d_strata_sample(y_pos, grid_size_, y_strata_)};
    float dist_sq = get_nearest_neighbor_dist_sq(cand_sample.x,
                                                 cand_sample.y);
    if (dist_sq > max_dist_sq) {
      max_sample = cand_sample;
      max_dist_sq = dist_sq;
    }
  }
  set_sample(sample_index, max_sample);
}

void SampleSet::set_sample(const int i,
                           const Sample& sample) {
  (*samples_)[i] = sample;

  x_strata_[sample.x * n_] = true;
  y_strata_[sample.y * n_] = true;

  int x_pos = sample.x * dim_;
  int y_pos = sample.y * dim_;
  sample_grid_[y_pos*dim_ + x_pos] = &(*samples_)[i];
}

float dist_sq(float x1, float y1, float x2, float y2) {
  float x_diff = x2-x1;
  float y_diff = y2-y1;
  return (x_diff*x_diff)+(y_diff*y_diff);
}

float SampleSet::get_nearest_neighbor_dist_sq(const float x, const float y) {
  // This function works by using the sample grid, since we know that the points
  // are well-distributed with at most one point in each cell. Much easier than
  // using any of the other data structures!
  //
  // Anyway start with the cells that are adjacent to our current cell and each
  // loop iteration we move outwards. We keep a track of the "grid radius",
  // which is the radius of the circle contained within our squares. If the
  // nearest point falls within this radius, we know that the next outward shift
  // can't find any nearer points.
  const int x_pos = x * dim_;
  const int y_pos = y * dim_;
  float min_dist_sq = 1.0;
  for (int i = 1; i <= dim_; i++) {
    float grid_radius = grid_size_ * i;
    float grid_radius_sq = grid_radius * grid_radius;

    const int x_min = std::max(x_pos - i, 0);
    const int x_max = std::min(x_pos + i, dim_-1);
    const int y_min = std::max(y_pos - i, 0);
    const int y_max = std::min(y_pos + i, dim_-1);
    // Traverse top and bottom boundaries, including corners.
    for (int x_offset = x_min; x_offset <= x_max; x_offset++) {
      const Sample* top_pt = sample_grid_[y_max*dim_+x_offset];
      const Sample* bottom_pt = sample_grid_[y_min*dim_+x_offset];
      if (top_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(top_pt->x, top_pt->y, x, y));
      if (bottom_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(bottom_pt->x, bottom_pt->y, x, y));
    }
    // Traverse left and right sides, excluding corners (hence the +1, -1).
    for (int y_offset = y_min+1; y_offset <= y_max-1; y_offset++) {
      const Sample* left_pt = sample_grid_[y_offset*dim_+x_min];
      const Sample* right_pt = sample_grid_[y_offset*dim_+x_max];
      if (left_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(left_pt->x, left_pt->y, x, y));
      if (right_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(right_pt->x, right_pt->y, x, y));
    }

    if (min_dist_sq < grid_radius_sq) {
      break;
    }
  }

  return min_dist_sq;
}

void pick_subquadrant_with_balance(const int x_balance,
                                   const int y_balance,
                                   int* x_pos,
                                   int* y_pos) {
  // The idea with the balance is, if we've mostly picked left cells within
  // subquadrants, we want to pick a right cell, and vice versa, and the same
  // for top and bottom. If we've netted 2 left cells and one top cell, we
  // care more about picking a right cell than a bottom cell.
  //
  // In Christensen et. al, they precompute these choices for the grid, but I
  // *THINK* this greedy method works just as well? And IMO it's easier to
  // understand.
  if (x_balance == 0 && y_balance == 0) {
    // If both are zero, we want to pick truly randomly.
    if (uniform_rand() < 0.5) *x_pos = *x_pos ^ 1;
    else
      *y_pos = *y_pos ^ 1;
  } else {
    bool use_x_balance = abs(x_balance) > abs(y_balance);
    if (abs(x_balance) == abs(y_balance)) {
      // If they're equal, we randomly pick one to balance. This is better
      // than randomly moving in one direction, because sometimes both
      // balances equally dictate the same move.
      use_x_balance = uniform_rand() < 0.5;
    }
    if (use_x_balance) {
      bool balance_to_right = x_balance < 0;
      if ((*x_pos & 1) != balance_to_right) *x_pos = *x_pos ^ 1;
      else
        *y_pos = *y_pos ^ 1;
    } else {
      bool balance_to_up = y_balance < 0;
      if ((*y_pos & 1) != balance_to_up) *y_pos = *y_pos ^ 1;
      else
        *x_pos = *x_pos ^ 1;
    }
  }
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

    pick_subquadrant_with_balance(x_balance, y_balance, &x_pos, &y_pos);

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
