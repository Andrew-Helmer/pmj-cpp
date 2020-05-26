// Copyright 2020 Andrew Helmer
#include "sample_generation/pmj.h"

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

class SampleSet {
 public:
  explicit SampleSet(const int num_samples,
                     const int num_candidates)
                     : num_samples(num_samples),
                       num_candidates_(num_candidates) {
    samples_ = std::make_unique<Point[]>(num_samples);

    x_strata_.resize(num_samples);
    y_strata_.resize(num_samples);

    int grid_memory_size = 1;
    while (grid_memory_size < num_samples)
      grid_memory_size <<= 2;
    sample_grid_ = std::make_unique<const Point*[]>(grid_memory_size);
  }

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

  // Gets the squared distance of the nearest neighbor to the given point.
  double GetNearestNeighborDistSq(const Point& sample) const;

  std::unique_ptr<Point[]> samples_;

  // Vector bool is usually implemented as a bitset!
  std::vector<bool> x_strata_ {false};
  std::vector<bool> y_strata_ {false};

  // The sample grid is used for nearest neighbor lookups.
  std::unique_ptr<const Point*[]> sample_grid_;

  int n_ = 1;  // Number of samples in the next pass.
  bool is_power_of_4_ = true;  // Whether n is a power of 4.
  int dim_ = 1;  // Number of cells in one dimension in next pass, i.e. sqrt(n).
  double grid_size_ = 1.0;  // 1.0 / dim_

  // Number of candidates to use for best-candidate sampling.
  const int num_candidates_;
};

void SampleSet::SubdivideStrata() {
  const int old_n = n_;

  n_ *= 2;
  is_power_of_4_ = !is_power_of_4_;
  if (!is_power_of_4_) {
    dim_ *= 2;
    grid_size_ *= 0.5;
  }

  std::fill_n(sample_grid_.get(), n_, nullptr);
  std::fill_n(x_strata_.begin(), n_, 0);
  std::fill_n(y_strata_.begin(), n_, 0);
  for (int i = 0; i < old_n; i++) {
    const auto& sample = samples_[i];

    x_strata_[sample.x * n_] = true;
    y_strata_[sample.y * n_] = true;

    const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
    sample_grid_[y_pos*dim_ + x_pos] = &sample;
  }
}

// This generates a sample within the grid position, verifying that it doesn't
// overlap strata with any other sample.
double get_1d_strata_sample(const int pos,
                            const int n,
                           const double grid_size,
                           const std::vector<bool>& strata) {
  while (true) {
    double val = UniformRand(pos*grid_size, (pos+1)*grid_size);
    int strata_pos = val * n;
    if (!strata[strata_pos]) {
      return val;
    }
  }
}

void SampleSet::GenerateNewSample(const int sample_index,
                                  const int x_pos,
                                  const int y_pos) {
  Point best_candidate;
  double max_dist_sq = 0.0;
  for (int i = 0; i < num_candidates_; i++) {
    Point cand_sample =
        {get_1d_strata_sample(x_pos, n_, grid_size_, x_strata_),
         get_1d_strata_sample(y_pos, n_, grid_size_, y_strata_)};
    if (num_candidates_ > 1) {
      double dist_sq = GetNearestNeighborDistSq(cand_sample);
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

void SampleSet::AddSample(const int i,
                          const Point& sample) {
  samples_[i] = sample;

  x_strata_[sample.x * n_] = true;
  y_strata_[sample.y * n_] = true;

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

/*
 * This function chooses the subquads to use in between even and odd powers
 * of two. In Christensen et al., they use ox-plowing, but only balance in the
 * Y dimension. This uses a Hilbert curve and a couple of other heuristics to
 * balance both dimensions.
 */
void HilbertRotate(int n, int *x, int *y, int rx, int ry) {
  // This function was copied from Wikipedia.
  if (ry == 0) {
      if (rx == 1) {
          *x = n-1 - *x;
          *y = n-1 - *y;
      }

      // Swap x and y.
      int t  = *x;
      *x = *y;
      *y = t;
  }
}
void HilbertIndexToPos(int n, int d, int *x, int *y) {
  // This function was copied from Wikipedia.
  int rx, ry, s, t = d;
  *x = *y = 0;
  for (s = 1; s < n; s *= 2) {
      rx = 1 & (t/2);
      ry = 1 & (t ^ rx);
      HilbertRotate(s, x, y, rx, ry);
      *x += s * rx;
      *y += s * ry;
      t /= 4;
  }
}
std::vector<std::pair<int, int>> GetBalancedChoices(const SampleSet& sample_set,
                                                    const int n) {
  // We'll return out choices at the end.
  std::vector<std::pair<int, int>> choices(n);

  // Sample Set dimension is the sub_quadrant dimension, we want the quadrant
  // dimension.
  const int quad_dim = sample_set.dim() / 2;

  // First we want to get the subquadrant positions, and also the sampling order
  // from
  std::vector<int> first_cells(n*2);
  std::vector<int> quadrant_order(n);
  for (int i = 0; i < n; i++) {
    const auto& sample = sample_set.sample(i);
    int x_pos = sample.x * sample_set.dim();
    int y_pos = sample.y * sample_set.dim();
    const int quadrant_index = (y_pos / 2)*(quad_dim) + (x_pos / 2);
    first_cells[2*quadrant_index] = x_pos;
    first_cells[2*quadrant_index+1] = y_pos;
    quadrant_order[quadrant_index] = i;
  }

  std::vector<int> choice_balance_x(quad_dim);
  std::vector<int> choice_balance_y(quad_dim);
  std::vector<int> num_visited_x(quad_dim);
  std::vector<int> num_visited_y(quad_dim);

  // This method doesn't alway work on the first try, though it usually does.
  // It pretty much never takes more than 3 or 4 tries, at it's pretty fast.
  for (int attempt = 0; attempt < 10; attempt++) {
    std::fill_n(choice_balance_x.begin(), quad_dim, 0);
    std::fill_n(choice_balance_y.begin(), quad_dim, 0);
    std::fill_n(num_visited_x.begin(), quad_dim, 0);
    std::fill_n(num_visited_y.begin(), quad_dim, 0);
    for (int i = 0; i < n; i++) {
      int col, row;
      HilbertIndexToPos(n, i, &col, &row);

      int quadrant_index = row*quad_dim + col;

      // The subquadrant positions of the first sample in this quadrant.
      int x_pos = first_cells[2*quadrant_index];
      int y_pos = first_cells[2*quadrant_index+1];

      // We'll either swap Y or X.
      bool swap_x = false;

      int balance_x = choice_balance_x[col];
      int balance_y = choice_balance_y[row];
      int x_remaining = quad_dim - num_visited_x[col] - abs(balance_x);
      int y_remaining = quad_dim - num_visited_y[row] - abs(balance_y);

      if (x_remaining == y_remaining && abs(balance_x) == abs(balance_y) &&
          balance_x != 0) {
        // If there's equal imbalance, we pick one of the balances randomly.
        if (UniformRand() < 0.5) balance_x = 0;
        else
          balance_y = 0;
      }

      if (x_remaining > y_remaining && balance_y != 0) {
        swap_x = (balance_y > 0) != (y_pos & 1);
      } else if (y_remaining > x_remaining && balance_x != 0) {
        swap_x = (balance_x > 0) == (x_pos & 1);
      } else if (abs(balance_x) > abs(balance_y)) {
        swap_x = (balance_x > 0) == (x_pos & 1);
      } else if (abs(balance_y) > abs(balance_x)) {
        swap_x = (balance_y > 0) != (y_pos & 1);
      } else {
        // balance_x == 0 && balance_y == 0 && x_remaining == y_remaining
        if (UniformRand() < 0.5) swap_x = true;
      }

      x_pos = swap_x ? (x_pos ^ 1) : x_pos;
      y_pos = !swap_x ? (y_pos ^ 1) : y_pos;

      choices[quadrant_order[quadrant_index]].first = x_pos;
      choices[quadrant_order[quadrant_index]].second = y_pos;

      num_visited_x[col]++;
      num_visited_y[row]++;

      choice_balance_x[col] += (x_pos & 1) ? 1 : -1;
      choice_balance_y[row] += (y_pos & 1) ? 1 : -1;
    }
    bool attempt_successful = true;
    for (int i = 0; i < quad_dim; i++) {
      if (choice_balance_x[i] != 0 || choice_balance_y[i] != 0) {
        attempt_successful = false;
      }
    }
    if (n == 1 || attempt_successful) {
      break;
    }
    std::cout << "Balance unsuccessful \n";
  }

  return choices;
}

std::unique_ptr<Point[]> GenerateSamples(
    const int num_samples,
    const int num_candidates) {
  SampleSet sample_set(num_samples, num_candidates);

  // Generate first sample.
  sample_set.GenerateNewSample(0, 0, 0);

  int n = 1;
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

    // Now we generate samples in the remaining subquadrants.
    sample_set.SubdivideStrata();

    // We want to make balanced choices here regarding which subquadrants to
    // use, so we precompute them in a special way (see above).
    auto sub_quad_choices = GetBalancedChoices(sample_set, n);
    for (int i = 0; i < n && 2*n+i < num_samples; i++) {
      sample_set.GenerateNewSample(
          2*n+i, sub_quad_choices[i].first, sub_quad_choices[i].second);
    }

    for (int i = 0; i < n && 3*n+i < num_samples; i++) {
      // Get the one diagonally opposite to the one we just got.
      sample_set.GenerateNewSample(
          3*n+i, sub_quad_choices[i].first ^ 1, sub_quad_choices[i].second ^ 1);
    }

    n *= 4;
  }

  return sample_set.ReleaseSamples();
}

}  // namespace

std::unique_ptr<Point[]> GetProgMultiJitteredSamples(
    const int num_samples) {
  return GenerateSamples(num_samples, 1);
}

std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoise(
    const int num_samples) {
  return GenerateSamples(num_samples, 10);
}

}  // namespace pmj
