// Copyright 2020 Andrew Helmer
#include "pmj.h"

#include <algorithm>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

namespace {

std::default_random_engine& get_rand_gen() {
  thread_local static std::random_device r;
  thread_local static std::default_random_engine gen(r());

  return gen;
}

float uniform_rand(float min, float max) {
  thread_local static std::uniform_real_distribution<float> uniform;

  std::uniform_real_distribution<float>::param_type param(min, max);

  return uniform(get_rand_gen(), param);
}
float uniform_rand() {
  return uniform_rand(0.0, 1.0);
}

std::pair<float, float> random_sample(
    float min_x, float max_x, float min_y, float max_y) {
  return {uniform_rand(min_x, max_x), uniform_rand(min_y, max_y)};
}

void prog_jittered_samples_quadrant(
    const std::pair<float, float>& sample,
    float min_x,
    float max_x,
    float min_y,
    float max_y,
    std::vector<std::pair<float, float>>* samples) {
  const float mid_x = 0.5*(min_x + max_x);
  const float mid_y = 0.5*(min_y + max_y);

  // Get quadrant of sample.
  bool is_left = sample.first < mid_x;
  bool is_bottom = sample.second < mid_y;

  // Generate diagonally opposite.
  samples->push_back(random_sample(is_left ? mid_x : min_x,
                                   is_left ? max_x : mid_x,
                                   is_bottom ? mid_y : min_y,
                                   is_bottom ? max_y : mid_y));
  // Switch value.
  if (uniform_rand() < 0.5) {
    is_left = !is_left;
  } else {
    is_bottom = !is_bottom;
  }

  samples->push_back(random_sample(is_left ? mid_x : min_x,
                                   is_left ? max_x : mid_x,
                                   is_bottom ? mid_y : min_y,
                                   is_bottom ? max_y : mid_y));

  // Switch both to do the diagonal of the last one.
  is_bottom = !is_bottom;
  is_left = !is_left;

  samples->push_back(random_sample(is_left ? mid_x : min_x,
                                   is_left ? max_x : mid_x,
                                   is_bottom ? mid_y : min_y,
                                   is_bottom ? max_y : mid_y));
}
}  // namespace

std::vector<std::pair<float, float>> prog_jittered_samples(
    const int num_samples) {
  std::vector<std::pair<float, float>> samples;

  // Generate first sample randomly
  samples.push_back(random_sample(0, 1, 0, 1));

  float grid_size = 1.0;
  while (samples.size() < num_samples) {
    std::vector<int> sample_indices(samples.size());
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    std::shuffle(sample_indices.begin(), sample_indices.end(), get_rand_gen());

    for (const int sample_index : sample_indices) {
      const auto& sample = samples[sample_index];

      float min_x = trunc(sample.first / grid_size)*grid_size;
      float min_y = trunc(sample.second / grid_size)*grid_size;

      prog_jittered_samples_quadrant(
          sample, min_x, min_x + grid_size, min_y, min_y + grid_size, &samples);

      if (samples.size() >= num_samples) {
        // At most we've done 3 too many.
        samples.resize(num_samples);
        break;
      }
    }
    grid_size *= 0.5;
  }

  return samples;
}

std::vector<std::pair<float, float>> prog_mj_samples(
    const int num_samples) {
  std::vector<std::pair<float, float>> samples = {
      {0.75, 0.75}, {0.25, 0.25}, {0.75, 0.25}, {0.25, 0.75}};
  return samples;
}

std::vector<std::pair<float, float>> prog_mj_samples_blue_noise(
    const int num_samples) {
  std::vector<std::pair<float, float>> samples = {
      {0.75, 0.75}, {0.25, 0.25}, {0.75, 0.25}, {0.25, 0.75}};
  return samples;
}
