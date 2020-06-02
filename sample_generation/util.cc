/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This file implements a few functions that are useful for multiple sampling
 * sequences, especially random number generation and nearest neighbor search
 * on a grid.
 */
#include "sample_generation/util.h"

#include <algorithm>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#include "sample_generation/pj.h"
#include "sample_generation/pmj.h"
#include "sample_generation/pmj02.h"
#include "sample_generation/pmj02bn.h"

namespace pmj {

thread_local static std::random_device r;
thread_local static std::default_random_engine gen(r());

double UniformRand(double min, double max) {
  thread_local static std::uniform_real_distribution<double> uniform;

  std::uniform_real_distribution<double>::param_type param(min, max);

  return uniform(gen, param);
}

int UniformInt(int min, int max) {
  thread_local static std::uniform_int_distribution<int> uniform;

  std::uniform_int_distribution<int>::param_type param(min, max);

  return uniform(gen, param);
}

namespace {
inline double GetToroidalDistSq(double x1, double y1, double x2, double y2) {
  double x_diff = abs(x2-x1);
  if (x_diff > 0.5) x_diff = 1.0 - x_diff;
  double y_diff = abs(y2-y1);
  if (y_diff > 0.5) y_diff = 1.0 - y_diff;

  return (x_diff*x_diff)+(y_diff*y_diff);
}

inline int WrapIndex(const int index,
                     const int limit) {
  if (index < 0) return index+limit;
  if (index >= limit) return index-limit;
  return index;
}

inline void UpdateMinDistSq(
    const Point& candidate, const Point& point, double* min_dist_sq) {
  double dist_sq =
      GetToroidalDistSq(point.x, point.y, candidate.x, candidate.y);
  if (dist_sq < *min_dist_sq) {
    *min_dist_sq = dist_sq;
  }
}

double GetNearestNeighborDistSq(const Point& sample,
                                const Point* sample_grid[],
                                const int dim,
                                const double max_min_dist_sq) {
  // This function works by using the sample grid, since we know that the points
  // are well-distributed with at most one point in each cell.
  //
  // Start with the cells that are adjacent to our current cell and each
  // loop iteration we move outwards. We keep a track of the "grid radius",
  // which is the radius of the circle contained within our squares. If the
  // nearest point falls within this radius, we know that the next outward shift
  // can't find any nearer points.
  //
  // We do wrap around cells, and compute toroidal distances.
  const int x_pos = sample.x * dim;
  const int y_pos = sample.y * dim;

  double min_dist_sq = 2.0;
  const double grid_size = 1.0 / dim;
  for (int i = 1; i <= dim; i++) {
    // We add sqrt(0.5) to account for the fact that the point might not be
    // in the center of its own square.
    double grid_radius = grid_size * (i + 0.7072);
    double grid_radius_sq = grid_radius * grid_radius;

    const int x_min = x_pos - i;
    const int x_max = x_pos + i;
    const int y_min = y_pos - i;
    const int y_max = y_pos + i;
    // Traverse top and bottom boundaries, including corners.
    for (int x_offset = x_min; x_offset <= x_max; x_offset++) {
      int wrapped_x_offset = WrapIndex(x_offset, dim);
      int wrapped_y_min = WrapIndex(y_min, dim);
      int wrapped_y_max = WrapIndex(y_max, dim);

      const Point* top_pt = sample_grid[wrapped_y_max*dim + wrapped_x_offset];
      const Point* bottom_pt =
          sample_grid[wrapped_y_min*dim + wrapped_x_offset];
      if (top_pt != nullptr) {
        min_dist_sq = std::min(min_dist_sq,
                               GetToroidalDistSq(top_pt->x, top_pt->y,
                                           sample.x, sample.y));
        UpdateMinDistSq(sample, *top_pt, &min_dist_sq);
      }
      if (bottom_pt != nullptr) {
        UpdateMinDistSq(sample, *bottom_pt, &min_dist_sq);
      }
    }
    // Traverse left and right sides, excluding corners (hence the +1, -1).
    for (int y_offset = y_min+1; y_offset <= y_max-1; y_offset++) {
      int wrapped_y_offset = WrapIndex(y_offset, dim);
      int wrapped_x_min = WrapIndex(x_min, dim);
      int wrapped_x_max = WrapIndex(x_max, dim);

      const Point* left_pt = sample_grid[wrapped_y_offset*dim + wrapped_x_min];
      const Point* right_pt = sample_grid[wrapped_y_offset*dim + wrapped_x_max];
      if (left_pt != nullptr) {
        UpdateMinDistSq(sample, *left_pt, &min_dist_sq);
      }
      if (right_pt != nullptr) {
        UpdateMinDistSq(sample, *right_pt, &min_dist_sq);
      }
    }

    if (min_dist_sq < grid_radius_sq ||
        min_dist_sq < max_min_dist_sq) {
      break;
    }
  }

  return min_dist_sq;
}
}  // namespace

Point GetBestCandidateOfSamples(const std::vector<Point>& candidates,
                                const Point* sample_grid[],
                                const int dim) {
  // Hypothetically, it could be faster to search all the points in parallel,
  // culling points as we go, but a naive implementation of this was only a tiny
  // bit faster, and the code was uglier, so we'll leave it for now.
  Point best_candidate;
  double max_min_dist_sq = 0.0;

  for (int i = 0; i < candidates.size(); i++) {
    Point cand_sample = candidates[i];
    double dist_sq =
        GetNearestNeighborDistSq(cand_sample,
                                 sample_grid,
                                 dim,
                                 max_min_dist_sq);
    if (dist_sq > max_min_dist_sq) {
      best_candidate = cand_sample;
      max_min_dist_sq = dist_sq;
    }
  }

  if (max_min_dist_sq_out != nullptr) {
    *max_min_dist_sq_out = max_min_dist_sq;
  }

  return best_candidate;
}

sample_fn GetSamplingFunction(const std::string& algorithm) {
  static const std::unordered_map<std::string, sample_fn> kAlgorithmMap = {
    {"uniform", &GetUniformRandomSamples},
    {"pj", &GetProgJitteredSamples},
    {"pmj", &GetProgMultiJitteredSamples},
    {"pmjbn", &GetProgMultiJitteredSamplesWithBlueNoise},
    {"pmj02", &GetPMJ02Samples},
    {"pmj02bn", &GetPMJ02SamplesWithBlueNoise},
    {"pmj02bn-2", &GetPMJ02SamplesWithBlueNoiseAttempts},
    /* Experimental/Explicit Algorithms */
    {"pmj-random", &GetProgMultiJitteredSamplesRandom},
    {"pmj-oxplowing", &GetProgMultiJitteredSamplesOxPlowing},
    {"pmj02-oxplowing", &GetProgMultiJitteredSamplesOxPlowing},
    {"pmj02-no-balance", &GetPMJ02SamplesNoBalance},
  };

  auto find_iterator = kAlgorithmMap.find(algorithm);
  if (find_iterator == kAlgorithmMap.end())
    throw std::invalid_argument(algorithm + " is not a valid algorithm.");

  return find_iterator->second;
}

/*
 * This is kind of like a binary tree shuffle. Easy to think about for 4 points.
 * We can swap points 1 and 2, and we can swap 3 and 4, and we can swap the
 * pairs (1,2) and (3,4), but we can't do anything else. So 2143 is a possible
 * sequence, but 1324 is not. This works because in a PMJ(0,2) sequence, every
 * multiple of a power of two is also a valid PMJ(0,2) sequence, at least if
 * it's constructed using the ShuffleSwap subquadrant selection.
 */
std::vector<const Point*> ShufflePMJ02Sequence(const pmj::Point points[],
                                               const int n) {
  assert((n & (n - 1)) == 0);  // This function only works for powers of two.
  std::vector<const Point*> shuffled_points(n);
  int random_encode = UniformInt(0, n-1);
  for (int i = 0; i < n; i++) {
    shuffled_points[i] = &points[i^random_encode];
  }

  return shuffled_points;
}

std::unique_ptr<Point[]> GetUniformRandomSamples(
    const int num_samples) {
  auto samples = std::make_unique<Point[]>(num_samples);

  for (int i = 0; i < num_samples; i++) {
    samples[i] = {UniformRand(), UniformRand()};
  }

  return samples;
}

}  // namespace pmj

