// Copyright 2020 Andrew Helmer
#include <iostream>
#include <memory>
#include <random>
#include <signal.h>
#include <utility>
#include <vector>

#include "pj.h"
#include "pmj.h"

int main() {
  // raise(SIGSTOP);
  std::unique_ptr<std::vector<pmj::Point>> samples =
      pmj::get_best_candidate_pmj_samples(4096);

  std::cout << "[";
  for (const auto& sample : (*samples)) {
    std::cout << "(" << sample.x << ", " << sample.y << "),\n";
  }
  std::cout << "]";
}
