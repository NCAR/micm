// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "robertson_system.hpp"

#include <iostream>

int main()
{
  auto sys = robertson::MakeSystem();
  std::cout << "reactions: " << sys.processes.size()
            << " B0=" << robertson::ConsistentB(0.04, 3e7, 1e4, 1.0, 0.0) << "\n";
  return 0;
}
