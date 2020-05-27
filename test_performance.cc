// Copyright 2020 Andrew Helmer.
#include <iostream>
#include <chrono>
#include <string>
#include <vector>

#include "absl/strings/str_split.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "sample_generation/pj.h"
#include "sample_generation/pmj.h"
#include "sample_generation/pmj02.h"

ABSL_FLAG(int, num_samples, 1024,
    "The number of samples to generate in each run.");
ABSL_FLAG(int, runs, 8, "The number of runs to make for each algorithm");
ABSL_FLAG(std::string, algorithms, "all",
    "Comma-separated list of algorithms to run. Use 'all' for all. Options are"
    "'pj', 'pmj', 'pmjbn', 'pmj02', 'pmj02bn'");

namespace chrono = std::chrono;

using std::string;
using std::vector;

typedef std::unique_ptr<pmj::Point[]> (*sample_f)(int);

sample_f GetSamplingFunction(const string& algorithm) {
  return algorithm == "pj" ? &pmj::GetProgJitteredSamples :
      algorithm == "pmj" ? &pmj::GetProgMultiJitteredSamples :
      algorithm == "pmjbn" ?
          &pmj::GetProgMultiJitteredSamplesWithBlueNoise :
      algorithm == "pmj02" ? &pmj::GetPMJ02Samples :
      algorithm == "pmj02bn" ? &pmj::GetPMJ02SamplesWithBlueNoise
      : throw std::invalid_argument(algorithm + " is not a valid algorithm.");
}

int main(int argc, char *argv[]) {
  absl::ParseCommandLine(argc, argv);

  vector<string> algorithms_list = {"pmj02"};
  if (absl::GetFlag(FLAGS_algorithms) == "all") {
    algorithms_list = {"pj", "pmj", "pmjbn", "pmj02", "pmj02bn"};
  } else {
    algorithms_list = absl::StrSplit(absl::GetFlag(FLAGS_algorithms), ',');
  }

  for (const string& algorithm : algorithms_list) {
    sample_f sample_func = GetSamplingFunction(algorithm);
    chrono::duration<double> total_time;
    for (int i = 0; i < absl::GetFlag(FLAGS_runs); i++) {
      auto start = chrono::high_resolution_clock::now();
      auto samples = (*sample_func)(absl::GetFlag(FLAGS_num_samples));
      auto end = chrono::high_resolution_clock::now();
      total_time += (end - start);
    }

    double samples_per_second = 1000000.0 * (absl::GetFlag(FLAGS_num_samples) * absl::GetFlag(FLAGS_runs)) /
        (chrono::duration_cast<chrono::microseconds>(total_time).count());

    std::cout << algorithm << ": " << samples_per_second << " samples/s\n";
  }

  return 0;
}
