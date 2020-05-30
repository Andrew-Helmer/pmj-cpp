/*
 * Copyright (C) Andrew Helmer 2020.
 *
 * Licensed under MIT Open-Source License: see LICENSE. If you use this code, or
 * you generate sample sets that you use, I'd appreciate a credit in the source
 * code of your software. Just my name and/or a link to the GitHub project.
 * Thanks!
 *
 * This file implements a few functions that are useful for multiple sampling
 * sequences, especially random number generation and nearest neighbor search
 * on a grid.
 */
#include "sample_generation/util.h"

#include <algorithm>
#include <random>
#include <utility>
#include <vector>

#include "sample_generation/pj.h"
#include "sample_generation/pmj.h"
#include "sample_generation/pmj02.h"

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
inline double DistSq(double x1, double y1, double x2, double y2) {
  double x_diff = x2-x1, y_diff = y2-y1;
  return (x_diff*x_diff)+(y_diff*y_diff);
}

double GetNearestNeighborDistSq(const Point& sample,
                                const Point* sample_grid[],
                                const int dim) {
  // This function works by using the sample grid, since we know that the points
  // are well-distributed with at most one point in each cell. Not the fastest
  // way to do this, but easy to implement and not horribly slow.
  //
  // Anyway start with the cells that are adjacent to our current cell and each
  // loop iteration we move outwards. We keep a track of the "grid radius",
  // which is the radius of the circle contained within our squares. If the
  // nearest point falls within this radius, we know that the next outward shift
  // can't find any nearer points.
  const int x_pos = sample.x * dim;
  const int y_pos = sample.y * dim;
  double min_dist_sq = sqrt(2.0);
  const double grid_size = 1.0 / dim;
  for (int i = 1; i <= dim; i++) {
    // We add sqrt(0.5) to account for the fact that the point might not be
    // in the center of its own square.
    double grid_radius = grid_size * (i + 0.7072);
    double grid_radius_sq = grid_radius * grid_radius;

    const int x_min = std::max(x_pos - i, 0);
    const int x_max = std::min(x_pos + i, dim-1);
    const int y_min = std::max(y_pos - i, 0);
    const int y_max = std::min(y_pos + i, dim-1);
    // Traverse top and bottom boundaries, including corners.
    for (int x_offset = x_min; x_offset <= x_max; x_offset++) {
      const Point* top_pt = sample_grid[y_max*dim+x_offset];
      const Point* bottom_pt = sample_grid[y_min*dim+x_offset];
      if (top_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 DistSq(top_pt->x, top_pt->y,
                                        sample.x, sample.y));
      if (bottom_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 DistSq(bottom_pt->x, bottom_pt->y,
                                        sample.x, sample.y));
    }
    // Traverse left and right sides, excluding corners (hence the +1, -1).
    for (int y_offset = y_min+1; y_offset <= y_max-1; y_offset++) {
      const Point* left_pt = sample_grid[y_offset*dim+x_min];
      const Point* right_pt = sample_grid[y_offset*dim+x_max];
      if (left_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 DistSq(left_pt->x, left_pt->y,
                                        sample.x, sample.y));
      if (right_pt != nullptr)
          min_dist_sq = std::min(min_dist_sq,
                                 DistSq(right_pt->x, right_pt->y,
                                        sample.x, sample.y));
    }

    if (min_dist_sq < grid_radius_sq) {
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
  double max_dist_sq = 0.0;

  for (int i = 0; i < candidates.size(); i++) {
    Point cand_sample = candidates[i];
    double dist_sq =
        GetNearestNeighborDistSq(cand_sample, sample_grid, dim);
    if (dist_sq > max_dist_sq) {
      best_candidate = cand_sample;
      max_dist_sq = dist_sq;
    }
  }

  return best_candidate;
}

sample_fn GetSamplingFunction(const std::string& algorithm) {
  return algorithm == "pj" ? &pmj::GetProgJitteredSamples :
      algorithm == "pmj" ? &pmj::GetProgMultiJitteredSamples :
            algorithm == "pmjbn" ?
          &pmj::GetProgMultiJitteredSamplesWithBlueNoise :
      algorithm == "pmj02" ? &pmj::GetPMJ02Samples :
      algorithm == "pmj02bn" ? &pmj::GetPMJ02SamplesWithBlueNoise :
      /* Experimental/Explicit Algorithms */
      algorithm == "pmj-oxplowing" ?
          &pmj::GetProgMultiJitteredSamplesOxPlowing :
      algorithm == "pmj-random" ? &pmj::GetProgMultiJitteredSamplesRandom :
      algorithm == "pmjbn-oxplowing" ?
          &pmj::GetProgMultiJitteredSamplesWithBlueNoiseOxPlowing :
      algorithm == "pmjbn-random" ?
          &pmj::GetProgMultiJitteredSamplesWithBlueNoiseRandom :
      algorithm == "pmj02-oxplowing" ? &pmj::GetPMJ02SamplesOxPlowing :
      algorithm == "pmj02-no-balance" ? &pmj::GetPMJ02SamplesNoBalance :
      throw std::invalid_argument(algorithm + " is not a valid algorithm.");
}

}  // namespace pmj

