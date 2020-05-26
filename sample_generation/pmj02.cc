// Copyright 2020 Andrew Helmer
#include "sample_generation/pmj02.h"

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

using std::vector;

class SampleSet {
 public:
  explicit SampleSet(const int num_samples,
                     const int num_candidates)
                     : num_samples(num_samples),
                       num_candidates_(num_candidates) {
    samples_ = std::make_unique<Point[]>(num_samples);

    int grid_memory_size = 1;
    while (grid_memory_size < num_samples)
      grid_memory_size <<= 2;
    sample_grid_.resize(grid_memory_size);
  }

  void PickSubquadrantWithBalance(const int sample_index,
                                     int* x_pos,
                                     int* y_pos);
  void GenerateFirstSample();
  // This generates a new sample at the current index, given the X position
  // and Y position of the subquadrant. It won't generate a new sample in an
  // existing strata.
  void GenerateNewSample(const int sample_index,
                         const int x_pos,
                         const int y_pos);

  // This function should be called after every power of 2 samples.
  void SubdivideStrata();

  // Get all the samples at the end.
  std::unique_ptr<Point[]> ReleaseSamples() {
    auto samples = std::make_unique<Point[]>(num_samples);
    samples.swap(samples_);
    return samples;
  }

  const Point& sample(const int sample_index) const {
    return samples_[sample_index];
  }
  const int dim() const { return dim_; }

  const int num_samples;

 private:
  // Adds a new point at index i. Updates the necessary data structures.
  void AddSample(const int i, const Point& sample);

  // Given a sample, sets all the correct strata to true.
  void UpdateStrata(const Point& sample, const int partial_strata_index = -1);

  Point GetCandidateSample(const int x_pos,
                           const int y_pos,
                           const int partial_strata_index = -1) const;

  bool IsStrataOccupied(const Point& sample,
                          const int partial_strata_index = -1) const;

  int GetPartialStrataIndex(const int sample_index) const;

  // Gets the squared distance of the nearest neighbor to the given point.
  double GetNearestNeighborDistSq(const Point& sample) const;

  std::unique_ptr<Point[]> samples_;

  // Contains all strata of elementary (0,2) intervals, except for the grid
  // where the stratum width and height are the same. Each value is true/false
  // for if a sample resides there.
  vector<vector<bool>> strata_ {};
  // This is the same as above, but it's only used for odd powers of two.
  vector<vector<bool>> partial_strata_[2] {};

  // The sample grid is used for nearest neighbor lookups.
  vector<const Point*> sample_grid_ {nullptr};

  int n_ = 1;  // Number of samples in the next pass.
  bool is_power_of_4_ = true;
  int dim_ = 1;  // Number of cells in one dimension in next pass, i.e. sqrt(n).
  double grid_size_ = 1.0;  // 1.0 / dim_

  // Number of candidates to use for best-candidate sampling.
  const int num_candidates_;
};

void SampleSet::SubdivideStrata() {
  const double old_n = n_;

  n_ *= 2;
  is_power_of_4_ = !is_power_of_4_;
  if (!is_power_of_4_) {
    dim_ *= 2;
    grid_size_ *= 0.5;
  }

  // For the first sample this is 1x1 (no strata). For sample 2 it's 1x2 and
  // 2x1. For samples 3-4 it's 4x1 and 1x4. For samples 5-8 it's 8x1, 2x4, 4x2,
  // 1x8. So it only goes up after reaching an odd power of two.
  if (!is_power_of_4_)
    strata_.resize(strata_.size()+2);

  std::fill(strata_.begin(), strata_.end(), vector<bool>(n_, false));

  // For 3-4 samples, the partial strata are both 1x1. For samples 9-12 and
  // 13-16, the partial strata are 4x1, 1x4. For samples 33-48 and 49-64, we
  // have 16x1, 8x2, 2x8, 1x16.
  if (is_power_of_4_ && n_ >= 16) {
    for (int i = 0; i < 2; i++) {
      partial_strata_[i].resize(partial_strata_[i].size()+2);
      std::fill(partial_strata_[i].begin(),
                partial_strata_[i].end(),
                vector<bool>(n_/4, false));
    }
  }

  std::fill_n(sample_grid_.begin(), old_n, nullptr);
  for (int i = 0; i < old_n; i++) {
    const auto& sample = samples_[i];

    UpdateStrata(sample);

    const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
    sample_grid_[y_pos*dim_ + x_pos] = &sample;
  }
}

// This generates a sample within the grid position, verifying that it doesn't
// overlap strata with any other sample.
Point SampleSet::GetCandidateSample(const int x_pos,
                                    const int y_pos,
                                    const int partial_strata_index) const {
  while (true) {
    // We first sample each dimension separately, and check only the 1D strata.
    // This significantly improves performance, with roughly a 2x speedup.
    Point sample;
    while (true) {
      sample.x = UniformRand(x_pos*grid_size_, (x_pos+1)*grid_size_);
      int strata_index = sample.x * n_;
      if (strata_[0][strata_index]) continue;
      break;
    }
    while (true) {
      sample.y = UniformRand(y_pos*grid_size_, (y_pos+1)*grid_size_);
      int strata_index = sample.y * n_;
      if (strata_.back()[strata_index]) continue;
      break;
    }

    if (!IsStrataOccupied(sample, partial_strata_index)) {
      return sample;
    }
  }
}

int SampleSet::GetPartialStrataIndex(const int sample_index) const {
  if (n_ <= 8) {
    // No partial strata yet.
    return -1;
  }
  // We use this to determine whether we're in an even or odd power of 2. Of
  // course there are other ways to do this.
  const int n_total = dim_*dim_;
  return sample_index < (n_total/2)
      ? -1
      : ((sample_index - n_total/2) < (n_total/4)
          ? 0
          : 1);
}

void SampleSet::GenerateFirstSample() {
  Point sample = {UniformRand(), UniformRand()};
  AddSample(0, sample);
}

void SampleSet::GenerateNewSample(const int sample_index,
                                  const int x_pos,
                                  const int y_pos) {
  Point best_candidate;
  double max_dist_sq = 0.0;

  for (int i = 0; i < num_candidates_; i++) {
    Point cand_sample = GetCandidateSample(
        x_pos, y_pos, GetPartialStrataIndex(sample_index));
    double dist_sq = GetNearestNeighborDistSq(cand_sample);
    if (dist_sq > max_dist_sq) {
      best_candidate = cand_sample;
      max_dist_sq = dist_sq;
    }
  }
  AddSample(sample_index, best_candidate);
}

void SampleSet::UpdateStrata(const Point& sample,
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

  // Only difference from above is n/4, and partial_strata_.
  if (partial_strata_index >= 0) {
    for (int i = 0, dim_x = n_ / 4, dim_y = 1;
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

bool SampleSet::IsStrataOccupied(const Point& sample,
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

  return false;

  // Only difference from above is n/4, and partial_strata_.
  if (partial_strata_index >= 0) {
      for (int i = 0, dim_x = n_ / 4, dim_y = 1;
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

void SampleSet::AddSample(const int i,
                          const Point& sample) {
  samples_[i] = sample;

  UpdateStrata(sample, GetPartialStrataIndex(i));

  const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
  sample_grid_[y_pos*dim_ + x_pos] = &(samples_[i]);
}

double dist_sq(double x1, double y1, double x2, double y2) {
  double x_diff = x2-x1, y_diff = y2-y1;
  return (x_diff*x_diff)+(y_diff*y_diff);
}

double SampleSet::GetNearestNeighborDistSq(const Point& sample) const {
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
  double min_dist_sq = 1.0;
  for (int i = 1; i <= dim_; i++) {
    double grid_radius = grid_size_ * i;
    double grid_radius_sq = grid_radius * grid_radius;

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
void SampleSet::PickSubquadrantWithBalance(const int sample_index,
                                           int* x_pos,
                                           int* y_pos) {
  const int partial_strata_index = GetPartialStrataIndex(sample_index);
  int x_pos1 = *x_pos ^ 1, y_pos1 = *y_pos;
  int x_pos2 = *x_pos , y_pos2 = *y_pos ^ 1;
  while (true) {
    Point sample1 = {UniformRand(x_pos1*grid_size_, (x_pos1+1)*grid_size_),
                     UniformRand(y_pos1*grid_size_, (y_pos1+1)*grid_size_)};
    Point sample2 = {UniformRand(x_pos2*grid_size_, (x_pos2+1)*grid_size_),
                     UniformRand(y_pos2*grid_size_, (y_pos2+1)*grid_size_)};
    if (!IsStrataOccupied(sample1, partial_strata_index)) {
      if (!IsStrataOccupied(sample2, partial_strata_index)) {
        bool use_subquad_1 = UniformRand() < 0.5;
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
    } else if (!IsStrataOccupied(sample2, partial_strata_index)) {
      *x_pos = x_pos2;
      *y_pos = y_pos2;
      return;
    }
  }
}

std::unique_ptr<Point[]> GenerateSamples(
    const int num_samples,
    const int num_candidates) {
  SampleSet sample_set(num_samples, num_candidates);

  // Generate first sample.
  sample_set.GenerateFirstSample();

  int n = 1;  // This n is not the same as the one in SampleSet!!
  while (n < num_samples) {
    sample_set.SubdivideStrata();

    // For every sample, we first generate the diagonally opposite one at the
    // current grid level.
    for (int i = 0; i < n && n+i < num_samples; i++) {
      const auto& sample = sample_set.sample(i);

      int x_pos = sample.x * sample_set.dim();
      int y_pos = sample.y * sample_set.dim();

      sample_set.GenerateNewSample(n+i, x_pos ^ 1, y_pos ^ 1);
      if (n+i >= num_samples) {
        break;
      }
    }

    sample_set.SubdivideStrata();

    // This represents the x_pos and y_pos of the subquadrants that we choose.
    int chosen_cells[n][2];
    for (int i = 0; i < n && 2*n+i < num_samples; i++) {
      const auto& sample = sample_set.sample(i);

      int x_pos = sample.x * sample_set.dim();
      int y_pos = sample.y * sample_set.dim();

      sample_set.PickSubquadrantWithBalance(2*n+i, &x_pos, &y_pos);

      sample_set.GenerateNewSample(2*n+i, x_pos, y_pos);

      chosen_cells[i][0] = x_pos;
      chosen_cells[i][1] = y_pos;
    }

    for (int i = 0; i < n && 3*n+i < num_samples; i++) {
      // Get the one diagonally opposite to the one we just got.
      sample_set.GenerateNewSample(
          3*n+i, chosen_cells[i][0] ^ 1, chosen_cells[i][1] ^ 1);
    }

    n *= 4;
  }

  return sample_set.ReleaseSamples();
}

}  // namespace

std::unique_ptr<Point[]> GetPMJ02Samples(
    const int num_samples) {
  return GenerateSamples(num_samples, 1);
}

std::unique_ptr<Point[]> GetPMJ02SamplesWithBlueNoise(
    const int num_samples) {
  return GenerateSamples(num_samples, 10);
}

}  // namespace pmj
