// Copyright 2020 Andrew Helmer.
#include <iostream>
#include <chrono>
#include <string>
#include <vector>

#include "absl/strings/str_split.h"
#include "gflags/gflags.h"
#include "sample_generation/pj.h"
#include "sample_generation/pmj.h"
#include "sample_generation/pmj02.h"

DEFINE_int32(n, 1024, "The number of samples to generate in each run.");
DEFINE_int32(runs, 8, "The number of runs to make for each algorithm");
DEFINE_string(algorithms, "pmj02",
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
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  vector<string> algorithms_list = {"pmj02"};
  if (FLAGS_algorithms == "all") {
    algorithms_list = {"pj, pmj, pmjbn, pmj02, pmj02bn"};
  } else {
    algorithms_list = absl::StrSplit(FLAGS_algorithms, ',');
  }

  for (const string& algorithm : algorithms_list) {
    sample_f sample_func = GetSamplingFunction(algorithm);
    chrono::duration<double> total_time;
    for (int i = 0; i < FLAGS_runs; i++) {
      auto start = chrono::high_resolution_clock::now();
      auto samples = (*sample_func)(FLAGS_n);
      auto end = chrono::high_resolution_clock::now();
      total_time += (end - start);
    }

    double samples_per_second = 1000000.0 * (FLAGS_n * FLAGS_runs) /
        (chrono::duration_cast<chrono::microseconds>(total_time).count());

    std::cout << algorithm << ": " << samples_per_second << " samples/s\n";
  }

  gflags::ShutDownCommandLineFlags();
  return 0;
}
