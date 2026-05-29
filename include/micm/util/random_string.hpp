// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <chrono>
#include <random>
#include <string>

namespace micm
{
  inline std::string GenerateRandomString()
  {
    auto now = std::chrono::system_clock::now();
    auto epoch = now.time_since_epoch();
    auto timestamp = std::chrono::duration_cast<std::chrono::seconds>(epoch).count();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 9999);
    int random_num = dist(gen);

    std::string random_num_str = std::to_string(random_num);
    std::string timestamp_str = std::to_string(timestamp);

    std::string result = timestamp_str + random_num_str;

    return result;
  }
}  // namespace micm