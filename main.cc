// Copyright 2020 Andrew Helmer
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include "pj.h"
#include "pmj.h"

int main() {
  std::vector<std::pair<float, float>> samples = prog_mj_samples(4096*16);

  for (const auto& sample : samples) {
    std::cout << "(" << sample.first << ", " << sample.second << ")\n";
  }
}
