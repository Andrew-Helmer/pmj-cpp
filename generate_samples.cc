/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This tool takes only two command-line arguments:
 * --algorithm can be one of "pj", "pmj", "pmjbn", "pmj02", or "pmj02bn". By
 * default it is pmj02.
 * --n is the number of samples, and it must be an integer greater than 0.
 *
 * Example usage with make:
 * $  make release
 * $  ./generate_samples --algorithm=pmj02bn --n=4096
 *
 * Example usage with bazel:
 * $  bazel run -c opt :generate_samples -- --algorithm=pmj02bn --n=4096
 */
#include <array>
#include <iostream>
#include <memory>
#include <random>
#include <utility>

#include "sample_generation/pj.h"
#include "sample_generation/pmj.h"
#include "sample_generation/pmj02.h"
#include "sample_generation/util.h"

using std::string;

namespace {
  void GetArguments(const int argc,
                    char* argv[],
                    int* n_samples,
                    std::string* algorithm) {
    string samples_str = "256";
    string algorithm_str = "pmj02";
    for (int i = 1; i < argc; ++i) {
      string arg = argv[i];
      if (arg == "--n") {
          if (i + 1 < argc) samples_str = argv[i+1];
      } else if (arg.rfind("--n=", 0) == 0) {
        samples_str = arg.substr(4);
      }

      if (arg == "--algorithm") {
          if (i + 1 < argc) algorithm_str = argv[i+1];
      } else if (arg.rfind("--algorithm=", 0) == 0) {
        algorithm_str = arg.substr(12);
      }
    }
    *n_samples = std::stoi(samples_str);
    if (*n_samples <= 0) {
      throw std::invalid_argument("--n must be positive.");
    }
    *algorithm = algorithm_str;
  }
}  // namespace

int main(int argc, char* argv[]) {
  int n_samples;
  string algorithm;
  GetArguments(argc, argv, &n_samples, &algorithm);
  auto* sample_func = pmj::GetSamplingFunction(algorithm);
  std::unique_ptr<pmj::Point[]> samples = (*sample_func)(n_samples);

  // Full double precision.
  std::cout.precision(17);

  for (int i = 0; i < n_samples; i++) {
    const auto& sample = samples[i];
    std::cout << "(" << sample.x << ", " << sample.y << "),\n";
  }
}
