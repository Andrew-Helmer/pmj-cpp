// Copyright 2020 Andrew Helmer
#include "sample_generation/pmj02.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <stack>
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
    sample_grid_ = std::make_unique<const Point*[]>(grid_memory_size);
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
    return std::move(samples_);
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
                           const vector<int>& valid_x_strata,
                           const vector<int>& valid_y_strata) const;

  bool IsStrataOccupied(const Point& sample,
                          const int partial_strata_index = -1) const;

  int GetPartialStrataIndex(const int sample_index) const;

  std::unique_ptr<Point[]> samples_;

  // Contains all strata of elementary (0,2) intervals. Each value is true/false
  // for if a sample resides there.
  vector<vector<bool>> strata_ {{false}};
  // This is the same as above, but it's only used for odd powers of two.
  vector<vector<bool>> partial_strata_[2] {};

  // The sample grid is used for nearest neighbor lookups.
  std::unique_ptr<const Point*[]> sample_grid_;

  int n_ = 1;  // Number of samples in the next pass.
  bool is_power_of_4_ = true;
  int dim_ = 1;  // Number of cells in one dimension in next pass, i.e. sqrt(n).

  // Number of candidates to use for best-candidate sampling.
  const int num_candidates_;
};

void SampleSet::SubdivideStrata() {
  const double old_n = n_;

  n_ *= 2;
  is_power_of_4_ = !is_power_of_4_;
  if (!is_power_of_4_) {
    dim_ *= 2;
  }

  // For the first sample this is 1x1. For sample 2 it's 1x2 and 2x1. For
  // samples 3-4 it's 4x1, 2x2, and 1x4. So every time it goes up by one.
  strata_.resize(strata_.size()+1);

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

  std::fill_n(sample_grid_.get(), old_n, nullptr);
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
                                    const vector<int>& valid_x_strata,
                                    const vector<int>& valid_y_strata) const {
  Point sample;

  int x_strata_index = valid_x_strata[UniformInt(0, valid_x_strata.size()-1)];
  int y_strata_index = valid_y_strata[UniformInt(0, valid_y_strata.size()-1)];

  double strata_width = 1.0 / n_;
  sample.x = UniformRand(strata_width*x_strata_index,
                         strata_width*(x_strata_index+1));
  sample.y = UniformRand(strata_width*y_strata_index,
                         strata_width*(y_strata_index+1));

  assert(sample.x >= 0.0 && sample.x < 1.0 && sample.y >= 0 && sample.y < 1.0);

  return sample;
}

int SampleSet::GetPartialStrataIndex(const int sample_index) const {
  return -1;

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

void GetXStrata(const int x_pos,
                const int y_pos,
                const int strata_index,
                const vector<vector<bool>>& strata,
                vector<int>* x_strata) {
  const int strata_x_dim = 1 << (strata.size() - strata_index - 1);
  const bool is_occupied =
      strata[strata_index][y_pos*strata_x_dim + x_pos];

  if (!is_occupied) {
    if (strata_index == 0) {
      // We're at the Nx1 leaf.
      x_strata->push_back(x_pos);
    } else {
      GetXStrata(x_pos * 2, y_pos / 2, strata_index - 1, strata, x_strata);
      GetXStrata(x_pos * 2 + 1, y_pos / 2, strata_index - 1, strata, x_strata);
    }
  }
}
void GetYStrata(const int x_pos,
                const int y_pos,
                const int strata_index,
                const vector<vector<bool>>& strata,
                vector<int>* y_strata) {
  const int strata_x_dim = 1 << (strata.size() - strata_index - 1);
  const bool is_occupied =
      strata[strata_index][y_pos*strata_x_dim + x_pos];

  if (!is_occupied) {
    if (strata_x_dim == 1) {
      // We're at the 1xN leaf.
      y_strata->push_back(y_pos);
    } else {
      GetYStrata(x_pos / 2, y_pos * 2, strata_index + 1, strata, y_strata);
      GetYStrata(x_pos / 2, y_pos * 2 + 1, strata_index + 1, strata, y_strata);
    }
  }
}

std::pair<vector<int>, vector<int>> GetValidStrata(
    const int x_pos, const int y_pos, const vector<vector<bool>>& strata) {
  std::pair<vector<int>, vector<int>> valid_strata = {{}, {}};

  if (strata.size() % 2 == 1) {
    GetXStrata(x_pos, y_pos, strata.size()/2, strata, &valid_strata.first);
    GetYStrata(x_pos, y_pos, strata.size()/2, strata, &valid_strata.second);
  } else {
    GetXStrata(x_pos, y_pos/2, strata.size()/2-1, strata, &valid_strata.first);
    GetYStrata(x_pos/2, y_pos, strata.size()/2, strata, &valid_strata.second);
  }

  return valid_strata;
}

void SampleSet::GenerateNewSample(const int sample_index,
                                  const int x_pos,
                                  const int y_pos) {
  Point best_candidate;
  double max_dist_sq = 0.0;

  int partial_strata_index = GetPartialStrataIndex(sample_index);

  const std::pair<vector<int>, vector<int>> valid_strata =
      GetValidStrata(x_pos, y_pos, strata_);

  for (int i = 0; i < num_candidates_; i++) {
    Point cand_sample = GetCandidateSample(
        x_pos, y_pos, valid_strata.first, valid_strata.second);
    if (num_candidates_ > 1) {
      double dist_sq =
          GetNearestNeighborDistSq(cand_sample, sample_grid_.get(), dim_);
      if (dist_sq > max_dist_sq) {
        best_candidate = cand_sample;
        max_dist_sq = dist_sq;
      }
    } else {
      best_candidate = cand_sample;
    }
  }
  AddSample(sample_index, best_candidate);
}

void SampleSet::UpdateStrata(const Point& sample,
                             const int partial_strata_index) {
  for (int i = 0, x_width = n_, y_width = 1;
       x_width >= 1;
       x_width /= 2, y_width *= 2, i++) {
    int x_pos = sample.x * x_width;
    int y_pos = sample.y * y_width;
    strata_[i][y_pos*x_width + x_pos] = true;
  }

  // Only difference from above is n/4, and partial_strata_.
  if (partial_strata_index >= 0) {
    for (int i = 0, x_width = n_ / 4, y_width = 1;
         x_width >= 1;
         x_width /= 2, y_width *= 2, i++) {
      int x_pos = sample.x * x_width;
      int y_pos = sample.y * y_width;

      partial_strata_[partial_strata_index][i][y_pos*x_width + x_pos] = true;
    }
  }
}

bool SampleSet::IsStrataOccupied(const Point& sample,
                                 const int partial_strata_index) const {
  for (int i = 0, x_width = n_, y_width = 1;
       x_width >= 1;
       x_width /= 2, y_width *= 2, i++) {
    int x_pos = sample.x * x_width;
    int y_pos = sample.y * y_width;

    if (strata_[i][y_pos*x_width + x_pos]) return true;
  }

  // Only difference from above is n/4, and partial_strata_.
  if (partial_strata_index >= 0) {
      for (int i = 0, x_width = n_ / 4, y_width = 1;
           x_width >= 1;
           x_width /= 2, y_width *= 2, i++) {
        if (x_width == y_width) {
          i--;
          continue;
        }
        int x_pos = sample.x * x_width;
        int y_pos = sample.y * y_width;

        if (partial_strata_[partial_strata_index][i][y_pos*x_width + x_pos])
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

// Modifies x_pos and y_pos to point to a new non-diagonal (adjacent)
// subquadrant, taking into account current balance.
void SampleSet::PickSubquadrantWithBalance(const int sample_index,
                                           int* x_pos,
                                           int* y_pos) {
  // No balance for now...
  bool use_subquad_1 = UniformRand() < 0.5;
  if (use_subquad_1) {
    *x_pos = *x_pos ^ 1;
    *y_pos = *y_pos;
  } else {
    *x_pos = *x_pos;
    *y_pos = *y_pos ^ 1;
  }
  return;

  // const int partial_strata_index = GetPartialStrataIndex(sample_index);
  // int x_pos1 = *x_pos ^ 1, y_pos1 = *y_pos;
  // int x_pos2 = *x_pos , y_pos2 = *y_pos ^ 1;
  // while (true) {
  //   Point sample1 = {UniformRand(x_pos1*grid_size_, (x_pos1+1)*grid_size_),
  //                    UniformRand(y_pos1*grid_size_, (y_pos1+1)*grid_size_)};
  //   Point sample2 = {UniformRand(x_pos2*grid_size_, (x_pos2+1)*grid_size_),
  //                    UniformRand(y_pos2*grid_size_, (y_pos2+1)*grid_size_)};
  //   if (!IsStrataOccupied(sample1, partial_strata_index)) {
  //     if (!IsStrataOccupied(sample2, partial_strata_index)) {
  //       bool use_subquad_1 = UniformRand() < 0.5;
  //       if (use_subquad_1) {
  //         *x_pos = x_pos1;
  //         *y_pos = y_pos1;
  //       } else {
  //         *x_pos = x_pos2;
  //         *y_pos = y_pos2;
  //       }
  //       return;
  //     }
  //     *x_pos = x_pos1;
  //     *y_pos = y_pos1;
  //     return;
  //   } else if (!IsStrataOccupied(sample2, partial_strata_index)) {
  //     *x_pos = x_pos2;
  //     *y_pos = y_pos2;
  //     return;
  //   }
  // }
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
