// Copyright 2020 Andrew Helmer
#ifndef UTIL_H_
#define UTIL_H_

#include <random>
#include <utility>
#include <vector>

std::default_random_engine& get_rand_gen();

float uniform_rand(float min, float max);

float uniform_rand();

std::pair<float, float> random_sample(
    float min_x, float max_x, float min_y, float max_y);

std::pair<float, float> get_sample(
    const int x_pos, const int y_pos, const float grid_size);

std::pair<float, float> get_diag_sample(
    const int x_pos, const int y_pos, const float grid_size);

#endif  // UTIL_H_
