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
  std::unique_ptr<std::vector<pmj::Sample>> samples =
      pmj::get_pmj_samples(256);

  std::cout << "[";
  for (const auto& sample : (*samples)) {
    std::cout << "(" << sample.x << ", " << sample.y << "),\n";
  }
  std::cout << "]";
}
