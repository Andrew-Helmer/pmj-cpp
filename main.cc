// Copyright 2020 Andrew Helmer
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include "pmj.h"

int main() {
  std::vector<std::pair<float, float>> samples = prog_jittered_samples(16);
  for (const auto& sample : samples) {
    std::cout << "(" << sample.first << ", " << sample.second << ")\n";
  }
}
