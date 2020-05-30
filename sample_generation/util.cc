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
inline double DistSqTorus(double x1, double y1, double x2, double y2) {
  double x_diff = abs(x2-x1);
  if (x_diff > 0.5) x_diff = 1.0 - x_diff;
  double y_diff = abs(y2-y1);
  if (y_diff > 0.5) y_diff = 1.0 - y_diff;

  return (x_diff*x_diff)+(y_diff*y_diff);
}

int Wrap(const int val,
         const int dim) {
  if (val < 0) return val+dim;
  if (val >= dim) return val-dim;
  return val;
}

double GetNearestNeighborDistSqTorus(const Point& sample,
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

    const int x_min = x_pos - i;
    const int x_max = x_pos + i;
    const int y_min = y_pos - i;
    const int y_max = y_pos + i;
    // Traverse top and bottom boundaries, including corners.
    for (int x_offset = x_min; x_offset <= x_max; x_offset++) {
      int wrapped_x_offset = Wrap(x_offset, dim);
      int wrapped_y_min = Wrap(y_min, dim);
      int wrapped_y_max = Wrap(y_max, dim);

      const Point* top_pt = sample_grid[wrapped_y_max*dim + wrapped_x_offset];
      const Point* bottom_pt =
          sample_grid[wrapped_y_min*dim + wrapped_x_offset];
      if (top_pt != nullptr) {
        min_dist_sq = std::min(min_dist_sq,
                               DistSqTorus(top_pt->x, top_pt->y,
                                           sample.x, sample.y));
      }
      if (bottom_pt != nullptr) {
        min_dist_sq = std::min(min_dist_sq,
                               DistSqTorus(bottom_pt->x, bottom_pt->y,
                                           sample.x, sample.y));
      }
    }
    // Traverse left and right sides, excluding corners (hence the +1, -1).
    for (int y_offset = y_min+1; y_offset <= y_max-1; y_offset++) {
      int wrapped_y_offset = Wrap(y_offset, dim);
      int wrapped_x_min = Wrap(x_min, dim);
      int wrapped_x_max = Wrap(x_max, dim);

      const Point* left_pt = sample_grid[wrapped_y_offset*dim + wrapped_x_min];
      const Point* right_pt = sample_grid[wrapped_y_offset*dim + wrapped_x_max];
      if (left_pt != nullptr) {
        min_dist_sq = std::min(min_dist_sq,
                               DistSqTorus(left_pt->x, left_pt->y,
                                           sample.x, sample.y));
      }
      if (right_pt != nullptr) {
        min_dist_sq = std::min(min_dist_sq,
                               DistSqTorus(right_pt->x, right_pt->y,
                                           sample.x, sample.y));
      }
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
        GetNearestNeighborDistSqTorus(cand_sample, sample_grid, dim);
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

