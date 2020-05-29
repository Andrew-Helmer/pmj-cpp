/*
 * Copyright (C) Andrew Helmer 2020.
 *
 * Licensed under MIT Open-Source License: see LICENSE.
 * If you use this code, or you generate samplesets that you use, I'd appreciate
 * a credit in the source code of your software. Just my name and/or a link to
 * the GitHub project. Thanks!
 *
 * Implementation of PMJ(0,2) sequences from
 * "Progressive Multi-Jittered Sample Sequences", Christensen et al. 2018, using
 * the algorithm from "Efficient Generation of Points that Satisfy
 * Two-Dimensional Elementary Intervals", Matt Pharr, 2019.
 *
 * This functionality is quite fast. On my 2017 Macbook Pro, 65536 samples takes
 * <50ms, i.e. it generates 1.43 million samples/sec, with -O3 compilation. Best
 * candidate sampling is slower at ~400,000 samples/sec with 20 candidates.
 * If you want to use the Best-Candidate samples in a production ray-tracer,
 * for example, probably better to precompute a bunch of tables and do lookups
 * into them. Also worth noting that the pmjbn algorithm (in pmj.cc) has much
 * better blue-noise characteristics than pmj02bn.
 */
#include "sample_generation/pmj02.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <stack>
#include <utility>
#include <vector>

#include "sample_generation/balance_util.h"
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

  void GenerateFirstSample();

  // This generates a new sample at the current index, given the X position
  // and Y position of the subquadrant. It won't generate a new sample in an
  // existing strata.
  void GenerateNewSample(const int sample_index,
                         const int x_pos,
                         const int y_pos);

  // This function should be called after every power of 2 samples.
  void SubdivideStrata();

  void ResetStrata(const int num_existing_samplesL);
  void InitSubsequenceStrata();

  // Get all the samples at the end.
  std::unique_ptr<Point[]> ReleaseSamples() {
    return std::move(samples_);
  }

  const Point& sample(const int sample_index) const {
    return samples_[sample_index];
  }
  const Point* samples() const {
    return samples_.get();
  }
  const int dim() const { return dim_; }

  const int num_samples;

 private:
  // Adds a new point at index i. Updates the necessary data structures.
  void AddSample(const int i, const Point& sample);

  // Given a sample, sets all the correct strata to true.
  void UpdateStrata(const Point& sample);

  Point GetCandidateSample(const int x_pos,
                           const int y_pos,
                           const vector<int>& valid_x_strata,
                           const vector<int>& valid_y_strata) const;

  std::unique_ptr<Point[]> samples_;

  // Contains all strata of elementary (0,2) intervals. Each value is true/false
  // for if a sample resides there.
  vector<vector<bool>> strata_ {{false}};

  vector<vector<bool>> subsequence_strata_ {{false}};

  // The sample grid is used for nearest neighbor lookups.
  std::unique_ptr<const Point*[]> sample_grid_;

  int n_ = 1;  // Number of samples in the next pass.
  bool is_power_of_4_ = true;
  bool use_subsequence_strata_ = false;
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

  ResetStrata(old_n);
}

void SampleSet::ResetStrata(const int num_existing_samples) {
  std::fill(strata_.begin(), strata_.end(), vector<bool>(n_, false));
  std::fill_n(sample_grid_.get(), n_, nullptr);

  use_subsequence_strata_ = false;

  for (int i = 0; i < num_existing_samples; i++) {
    const auto& sample = samples_[i];

    UpdateStrata(sample);

    const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
    sample_grid_[y_pos*dim_ + x_pos] = &sample;
  }
}

void SampleSet::InitSubsequenceStrata() {
  use_subsequence_strata_ = true;
  subsequence_strata_ = strata_;
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

void SampleSet::GenerateFirstSample() {
  Point sample = {UniformRand(), UniformRand()};
  AddSample(0, sample);
}

inline void GetXStrata(const int x_pos,
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
inline void GetYStrata(const int x_pos,
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

inline std::pair<vector<int>, vector<int>> GetValidStrata(
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

  vector<vector<bool>>* strata = &strata_;
  if (use_subsequence_strata_) {
    strata = &subsequence_strata_;
  }
  const std::pair<vector<int>, vector<int>>& valid_strata =
      GetValidStrata(x_pos, y_pos, *strata);

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

void SampleSet::UpdateStrata(const Point& sample) {
  vector<vector<bool>>* strata = &strata_;
  if (use_subsequence_strata_) {
    strata = &subsequence_strata_;
  }

  for (int i = 0, x_width = n_, y_width = 1;
       x_width >= 1;
       x_width /= 2, y_width *= 2, i++) {
    int x_pos = sample.x * x_width;
    int y_pos = sample.y * y_width;
    (*strata)[i][y_pos*x_width + x_pos] = true;

    if (use_subsequence_strata_ && n_ > 4 && x_width >= 2 && y_width >= 2) {
      (*strata)[i][(y_pos^1)*x_width + x_pos] = true;
      (*strata)[i][y_pos*x_width + (x_pos^1)] = true;
      (*strata)[i][(y_pos^1)*x_width + (x_pos^1)] = true;
    }
  }
}

void SampleSet::AddSample(const int i,
                          const Point& sample) {
  samples_[i] = sample;

  UpdateStrata(sample);

  const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
  sample_grid_[y_pos*dim_ + x_pos] = &(samples_[i]);
}

std::unique_ptr<Point[]> GenerateSamples(
    const int num_samples,
    const int num_candidates,
    const bool subsequence_stratification = true,
    const subquad_fn subquad_func = &GetSubQuadrantsConsistently) {
  SampleSet sample_set(num_samples, num_candidates);

  // A different balance function doesn't work with subsequence stratification,
  // it gets stuck.
  assert(subquad_func == (&GetSubQuadrantsConsistently) ||
         !subsequence_stratification);

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
    if (subsequence_stratification) {
      sample_set.InitSubsequenceStrata();
    }

    auto sub_quad_choices =
        (*subquad_func)(sample_set.samples(), sample_set.dim());
    for (int i = 0; i < n && 2*n+i < num_samples; i++) {
      sample_set.GenerateNewSample(
          2*n+i, sub_quad_choices[i].first, sub_quad_choices[i].second);
    }

    if (subsequence_stratification) {
      sample_set.ResetStrata(3*n);
      sample_set.InitSubsequenceStrata();
    }

    for (int i = 0; i < n && 3*n+i < num_samples; i++) {
      // Get the one diagonally opposite to the one we just got.
      sample_set.GenerateNewSample(
          3*n+i, sub_quad_choices[i].first ^ 1, sub_quad_choices[i].second^ 1);
    }

    n *= 4;
  }

  return sample_set.ReleaseSamples();
}

}  // namespace

std::unique_ptr<Point[]> GetPMJ02Samples(
    const int num_samples) {
  return GenerateSamples(num_samples, /*num_candidates=*/1);
}

std::unique_ptr<Point[]> GetPMJ02SamplesWithBlueNoise(
    const int num_samples) {
  return GenerateSamples(num_samples, /*num_candidates=*/20);
}

// These functions are just for experimentation, no reason to actually use them.
std::unique_ptr<Point[]> GetPMJ02SamplesNoBalance(
    const int num_samples) {
  return GenerateSamples(num_samples,
                         /*num_candidates=*/20,
                         /*subsequence_stratification=*/false,
                         /*subquad_func=*/&GetSubQuadrantsRandomly);
}
std::unique_ptr<Point[]> GetPMJ02SamplesOxPlowing(
    const int num_samples) {
  return GenerateSamples(num_samples,
                         /*num_candidates=*/20,
                         /*subsequence_stratification=*/false,
                         /*subquad_func=*/&GetSubQuadrantsOxPlowing);
}
}  // namespace pmj
