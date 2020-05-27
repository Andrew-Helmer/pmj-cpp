// Copyright 2020 Andrew Helmer
#include "sample_generation/pmj_util.h"

#include <utility>
#include <vector>

#include "sample_generation/util.h"

namespace pmj {
namespace {
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

}  // namespace

/*
 * This function chooses the subquads to use in between even and odd powers
 * of two. In Christensen et al., they use ox-plowing, but only balance in the
 * Y dimension. This uses a Hilbert curve and a couple of other heuristics to
 * balance both dimensions.
 */
std::vector<std::pair<int, int>> GetBalancedChoicesHilbert(
    const Point samples[],
    const int dim) {
  // Dim is the subquadrant dimension, we want the quadrant dimension.
  const int quad_dim = dim / 2;
  const int n = quad_dim*quad_dim;

  // We'll return out choices at the end.
  std::vector<std::pair<int, int>> choices(n);

  // First we want to get the subquadrant positions, and also the sampling order
  // from
  std::vector<int> first_cells(n*2);
  std::vector<int> quadrant_order(n);
  for (int i = 0; i < n; i++) {
    const auto& sample = samples[i];
    int x_pos = sample.x * dim;
    int y_pos = sample.y * dim;
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

      // We'll either swap X or Y.
      bool swap_x = false;

      int balance_x = choice_balance_x[col];
      int balance_y = choice_balance_y[row];

      const int x_surplus = (quad_dim - num_visited_x[col]) - abs(balance_x);
      const int y_surplus = (quad_dim - num_visited_y[row]) - abs(balance_y);

      if (x_surplus == y_surplus && abs(balance_x) == abs(balance_y) &&
          balance_x != 0) {
        // If there's equal imbalance, we pick one of the balances randomly.
        if (UniformRand() < 0.5) balance_x = 0;
        else
          balance_y = 0;
      }

      if (x_surplus > y_surplus && balance_y != 0) {
        swap_x = (balance_y > 0) != (y_pos & 1);
      } else if (y_surplus > x_surplus && balance_x != 0) {
        swap_x = (balance_x > 0) == (x_pos & 1);
      } else if (abs(balance_x) > abs(balance_y)) {
        swap_x = (balance_x > 0) == (x_pos & 1);
      } else if (abs(balance_y) > abs(balance_x)) {
        swap_x = (balance_y > 0) != (y_pos & 1);
      } else {
        // Everything is 0!
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
  }

  return choices;
}

}  // namespace pmj
