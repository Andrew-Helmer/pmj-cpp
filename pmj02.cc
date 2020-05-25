// Copyright 2020 Andrew Helmer
#include "pmj02.h"

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
  explicit SampleSet(const int num_samples,
                     const int num_candidates)
                     : num_samples(num_samples),
                       num_candidates_(num_candidates) {
    samples_ = std::make_unique<std::vector<Point>>();
    samples_->reserve(num_samples);
    samples_->push_back({0.0, 0.0});

    sample_grid_.reserve(num_samples);
  }

  void pick_subquadrant_with_balance(const int sample_index,
                                     int* x_pos,
                                     int* y_pos);

  // This generates a new sample at the current index, given the X position
  // and Y position of the subquadrant. It won't generate a new sample in an
  // existing strata.
  void create_new_sample(const int sample_index,
                         const int x_pos,
                         const int y_pos);

  // This function should be called after every power of 4 samples. E.g. after
  // 1 sample, 4, 16, etc. It divides all the strata for the next pass, and
  // figures out where the existing points lie.
  void subdivide_strata();

  void subdivide_strata_2x();

  // Get all the samples at the end.
  std::unique_ptr<std::vector<Point>> release_samples() {
    auto samples = std::make_unique<std::vector<Point>>();
    samples.swap(samples_);
    return samples;
  }

  const Point& sample(const int sample_index) const {
    return (*samples_)[sample_index];
  }
  const int dim() const { return dim_; }
  const float grid_size() const { return grid_size_; }

  const int num_samples;

 private:
  // Adds a new point at index i. Updates the necessary data structures.
  void set_sample(const int i, const Point& sample);

  // Given a sample, sets all the correct strata to true.
  void set_strata(const Point& sample, const int partial_strata_index = -1);

  Point get_candidate_sample(const int x_pos,
                             const int y_pos,
                             const int partial_strata_index = -1) const;

  bool is_strata_occupied(const Point& sample,
                          const int partial_strata_index = -1) const;

  int get_partial_strata_index(const int sample_index) const;

  // Gets the squared distance of the nearest neighbor to the given point.
  float get_nearest_neighbor_dist_sq(const Point& sample) const;

  std::unique_ptr<std::vector<Point>> samples_;

  // Contains all strata of elementary (0,2) intervals, except for the grid
  // where the stratum width and height are the same. Each value is true/false
  // for if a sample resides there.
  std::vector<std::vector<bool>> strata_ {};
  // This is the same as above, but it's only used for odd powers of two.
  std::vector<std::vector<bool>> partial_strata_[2] {};

  // The sample grid is used for nearest neighbor lookups.
  std::vector<const Point*> sample_grid_ {nullptr};

  int n_ = 1;  // Number of samples in the next pass.
  int dim_ = 1;  // Number of cells in one dimension in next pass, i.e. sqrt(n).
  float grid_size_ = 1.0;  // 1.0 / dim_

  // Number of candidates to use for best-candidate sampling.
  const int num_candidates_;
};

void SampleSet::subdivide_strata() {
  const float old_n = n_;

  n_ *= 2;
  dim_ *= 2;
  grid_size_ *= 0.5;

  samples_->resize(std::min(n_, num_samples));

  // For the first sample this is 1x1 (no strata). For sample 2 it's 1x2 and
  // 2x1. For samples 3-4 it's 4x1 and 1x4. For samples 5-8 it's 8x1, 2x4, 4x2,
  // 1x8. So it only goes up after reaching a power of 4.
  strata_.resize(strata_.size()+2);
  std::fill(strata_.begin(), strata_.end(), std::vector<bool>(n_, false));

  sample_grid_.resize(n_);
  std::fill(sample_grid_.begin(), sample_grid_.end(), nullptr);
  for (int i = 0; i < old_n; i++) {
    const auto& sample = (*samples_)[i];

    set_strata(sample);

    const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
    sample_grid_[y_pos*dim_ + x_pos] = &sample;
  }
}

void SampleSet::subdivide_strata_2x() {
  const int old_n = n_;
  n_ *= 2;

  samples_->resize(std::min(n_, num_samples));

  // The number of strata doesn't change, but the size of each strata does.
  // For instance, it might go from 1x2 to 2x4.
  std::fill(strata_.begin(), strata_.end(), std::vector<bool>(n_, false));

  // For 3-4 samples, the partial strata are 2x1 and 1x2. For 9-16, they are
  // 4x1 and 1x4. For 33-64, they are 16x1, 8x2, 2x8, 4x4.
  for (int i = 0; i < 2; i++) {
    partial_strata_[i].resize(partial_strata_[i].size()+2);
    std::fill(partial_strata_[i].begin(),
              partial_strata_[i].end(),
              std::vector<bool>(n_/2, false));
  }

  for (int i = 0; i < old_n; i++) {
    const auto& sample = (*samples_)[i];

    set_strata(sample);
  }
}

// This generates a sample within the grid position, verifying that it doesn't
// overlap strata with any other sample.
Point SampleSet::get_candidate_sample(const int x_pos,
                                      const int y_pos,
                                      const int partial_strata_index) const {
  while (true) {
    Point sample = {uniform_rand(x_pos*grid_size_, (x_pos+1)*grid_size_),
                    uniform_rand(y_pos*grid_size_, (y_pos+1)*grid_size_)};
    if (!is_strata_occupied(sample, partial_strata_index)) {
      return sample;
    }
  }
}

int SampleSet::get_partial_strata_index(const int sample_index) const {
  // We use this to determine whether we're in an even or odd power of 2. Of
  // course there are other ways to do this.
  const int n_total = dim_*dim_;
  return sample_index < (n_total/2)
      ? -1
      : ((sample_index - n_total/2) < (n_total/2)
          ? 0
          : 1);
}

void SampleSet::create_new_sample(const int sample_index,
                                  const int x_pos,
                                  const int y_pos) {
  Point best_candidate;
  float max_dist_sq = 0.0;

  for (int i = 0; i < num_candidates_; i++) {
    Point cand_sample = get_candidate_sample(
        x_pos, y_pos, get_partial_strata_index(sample_index));
    float dist_sq = get_nearest_neighbor_dist_sq(cand_sample);
    if (dist_sq > max_dist_sq) {
      best_candidate = cand_sample;
      max_dist_sq = dist_sq;
    }
  }
  set_sample(sample_index, best_candidate);
}

void SampleSet::set_strata(const Point& sample,
                           const int partial_strata_index) {
  for (int i = 0, dim_x = n_, dim_y = 1;
       dim_x >= 1;
       dim_x /= 2, dim_y *= 2, i++) {
    if (dim_x == dim_y) {
      i--;
      continue;
    }
    int x_pos = sample.x * dim_x;
    int y_pos = sample.y * dim_y;
    strata_[i][y_pos*dim_x + x_pos] = true;
  }

  // Only difference from above is n/2, and partial_strata_.
  if (partial_strata_index >= 0) {
    for (int i = 0, dim_x = n_ / 2, dim_y = 1;
         dim_x >= 1;
         dim_x /= 2, dim_y *= 2, i++) {
      if (dim_x == dim_y) {
        i--;
        continue;
      }
      int x_pos = sample.x * dim_x;
      int y_pos = sample.y * dim_y;

      partial_strata_[partial_strata_index][i][y_pos*dim_x + x_pos] = true;
    }
  }
}

bool SampleSet::is_strata_occupied(const Point& sample,
                                   const int partial_strata_index) const {
  for (int i = 0, dim_x = n_, dim_y = 1;
       dim_x >= 1;
       dim_x /= 2, dim_y *= 2, i++) {
    if (dim_x == dim_y) {
      i--;
      continue;
    }
    int x_pos = sample.x * dim_x;
    int y_pos = sample.y * dim_y;

    if (strata_[i][y_pos*dim_x + x_pos]) return true;
  }

  // Only difference from above is n/2, and partial_strata_.
  if (partial_strata_index >= 0) {
      for (int i = 0, dim_x = n_ / 2, dim_y = 1;
           dim_x >= 1;
           dim_x /= 2, dim_y *= 2, i++) {
        if (dim_x == dim_y) {
          i--;
          continue;
        }
        int x_pos = sample.x * dim_x;
        int y_pos = sample.y * dim_y;

        if (partial_strata_[partial_strata_index][i][y_pos*dim_x + x_pos])
          return true;
      }
  }
  return false;
}

void SampleSet::set_sample(const int i,
                           const Point& sample) {
  (*samples_)[i] = sample;

  set_strata(sample, get_partial_strata_index(i));

  const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
  sample_grid_[y_pos*dim_ + x_pos] = &(*samples_)[i];
}

float dist_sq(float x1, float y1, float x2, float y2) {
  float x_diff = x2-x1, y_diff = y2-y1;
  return (x_diff*x_diff)+(y_diff*y_diff);
}

float SampleSet::get_nearest_neighbor_dist_sq(const Point& sample) const {
  // This function works by using the sample grid, since we know that the points
  // are well-distributed with at most one point in each cell. Much easier than
  // using any of the other data structures!
  //
  // Anyway start with the cells that are adjacent to our current cell and each
  // loop iteration we move outwards. We keep a track of the "grid radius",
  // which is the radius of the circle contained within our squares. If the
  // nearest point falls within this radius, we know that the next outward shift
  // can't find any nearer points.
  const int x_pos = sample.x * dim_;
  const int y_pos = sample.y * dim_;
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
      const Point* top_pt = sample_grid_[y_max*dim_+x_offset];
      const Point* bottom_pt = sample_grid_[y_min*dim_+x_offset];
      if (top_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(top_pt->x, top_pt->y,
                                         sample.x, sample.y));
      if (bottom_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(bottom_pt->x, bottom_pt->y,
                                         sample.x, sample.y));
    }
    // Traverse left and right sides, excluding corners (hence the +1, -1).
    for (int y_offset = y_min+1; y_offset <= y_max-1; y_offset++) {
      const Point* left_pt = sample_grid_[y_offset*dim_+x_min];
      const Point* right_pt = sample_grid_[y_offset*dim_+x_max];
      if (left_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(left_pt->x, left_pt->y,
                                         sample.x, sample.y));
      if (right_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 dist_sq(right_pt->x, right_pt->y,
                                         sample.x, sample.y));
    }

    if (min_dist_sq < grid_radius_sq) {
      break;
    }
  }

  return min_dist_sq;
}

// Modifies x_pos and y_pos to point to a new non-diagonal (adjacent)
// subquadrant, taking into account current balance.
void SampleSet::pick_subquadrant_with_balance(const int sample_index,
                                              int* x_pos,
                                              int* y_pos) {
  const int partial_strata_index = get_partial_strata_index(sample_index);
  int x_pos1 = *x_pos ^ 1, y_pos1 = *y_pos;
  int x_pos2 = *x_pos , y_pos2 = *y_pos ^ 1;
  while (true) {
    Point sample1 = {uniform_rand(x_pos1*grid_size_, (x_pos1+1)*grid_size_),
                     uniform_rand(y_pos1*grid_size_, (y_pos1+1)*grid_size_)};
    Point sample2 = {uniform_rand(x_pos2*grid_size_, (x_pos2+1)*grid_size_),
                     uniform_rand(y_pos2*grid_size_, (y_pos2+1)*grid_size_)};
    if (!is_strata_occupied(sample1, partial_strata_index)) {
      if (!is_strata_occupied(sample2, partial_strata_index)) {
        // In this case we actually get candidate samples from both, and pick
        // the best based off that.
        float max_dist_sq = 0.0;

        bool use_subquad_1 = uniform_rand() < 0.5;
        for (int i = 0; i < num_candidates_; i++) {
          Point cand_sample1 = {uniform_rand(x_pos1*grid_size_, (x_pos1+1)*grid_size_),
                     uniform_rand(y_pos1*grid_size_, (y_pos1+1)*grid_size_)};
          Point cand_sample2 = {uniform_rand(x_pos2*grid_size_, (x_pos2+1)*grid_size_),
                           uniform_rand(y_pos2*grid_size_, (y_pos2+1)*grid_size_)};
          float dist_sq1 = get_nearest_neighbor_dist_sq(cand_sample1);
          float dist_sq2 = get_nearest_neighbor_dist_sq(cand_sample2);
          if (dist_sq1 > max_dist_sq) {
            max_dist_sq = dist_sq1;
            use_subquad_1 = true;
          }
          if (dist_sq2 > max_dist_sq) {
            max_dist_sq = dist_sq2;
            use_subquad_1 = false;
          }
        }

        if (use_subquad_1) {
          *x_pos = x_pos1;
          *y_pos = y_pos1;
        } else {
          *x_pos = x_pos2;
          *y_pos = y_pos2;
        }
        return;
      }
      *x_pos = x_pos1;
      *y_pos = y_pos1;
      return;
    } else if (!is_strata_occupied(sample2, partial_strata_index)) {
      *x_pos = x_pos2;
      *y_pos = y_pos2;
      return;
    }
  }
}

std::unique_ptr<std::vector<Point>> get_samples(
    const int num_samples,
    const int num_candidates) {
  SampleSet sample_set(num_samples, num_candidates);

  // Generate first sample.
  sample_set.create_new_sample(0, 0, 0);

  int n = 1;  // This n is not the same as the one in SampleSet!!
  while (n < num_samples) {
    sample_set.subdivide_strata();

    // For every sample, we first generate the diagonally opposite one at the
    // current grid level.
    for (int i = 0; i < n && n+i < num_samples; i++) {
      const auto& sample = sample_set.sample(i);

      int x_pos = sample.x * sample_set.dim();
      int y_pos = sample.y * sample_set.dim();

      sample_set.create_new_sample(n+i, x_pos ^ 1, y_pos ^ 1);
      if (n+i >= num_samples) {
        break;
      }
    }

    sample_set.subdivide_strata_2x();
    for (int i = 0; i < n && 2*n+i < num_samples; i++) {
      const auto& sample = sample_set.sample(i);

      int x_pos = sample.x * sample_set.dim();
      int y_pos = sample.y * sample_set.dim();

      sample_set.pick_subquadrant_with_balance(2*n+i, &x_pos, &y_pos);

      sample_set.create_new_sample(2*n+i, x_pos, y_pos);

      // Get the one diagonally opposite to the one we just got.
      if (3*n+i >= num_samples) continue;
      sample_set.create_new_sample(3*n+i, x_pos ^ 1, y_pos ^ 1);
    }

    n *= 4;
  }

  return sample_set.release_samples();
}

}  // namespace

std::unique_ptr<std::vector<Point>> get_pmj02_samples(
    const int num_samples) {
  return get_samples(num_samples, 1);
}

std::unique_ptr<std::vector<Point>> get_best_candidate_pmj02_samples(
    const int num_samples) {
  return get_samples(num_samples, 10);
}

}  // namespace pmj
