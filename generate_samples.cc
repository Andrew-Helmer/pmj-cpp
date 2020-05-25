// Copyright 2020 Andrew Helmer.
// This simple function takes only two command-line arguments, algorithm and
// n. Algorithm can be one of "pj", "pmj", "pmjbn", "pmj02", or "pmj02bn". By
// default it is pmj02. n must be an integer greater than zero.
// Example usage:
// $  make
// $  ./generate_samples --algorithm='pmj02bn' --n=4096
#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "pj.h"
#include "pmj.h"
#include "pmj02.h"

namespace {
  void GetArguments(const int argc,
                    char* argv[],
                    int* n_samples,
                    std::string* algorithm) {
    std::string samples_str = "256";
    std::string algorithm_str = "pmj02";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
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
  std::string algorithm;
  GetArguments(argc, argv, &n_samples, &algorithm);
  std::unique_ptr<std::vector<pmj::Point>> samples =
      algorithm == "pmj" ? pmj::get_pmj_samples(n_samples) :
      algorithm == "pmjbn" ? pmj::get_best_candidate_pmj_samples(n_samples) :
      algorithm == "pmj02" ? pmj::get_pmj02_samples(n_samples) :
      algorithm == "pmj02bn" ? pmj::get_best_candidate_pmj02_samples(n_samples)
      : throw std::invalid_argument(algorithm + " is not a valid algorithm.");

  std::cout << "[";
  for (const auto& sample : (*samples)) {
    std::cout << "(" << sample.x << ", " << sample.y << "),\n";
  }
  std::cout << "]";
}
