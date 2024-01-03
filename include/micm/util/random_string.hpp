// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <chrono>
#include <random>
#include <string>

namespace micm
{
  std::string generate_random_string()
  {
    auto now = std::chrono::system_clock::now();
    auto epoch = now.time_since_epoch();
    auto timestamp = std::chrono::duration_cast<std::chrono::seconds>(epoch).count();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 9999);
    int randomNum = dist(gen);

    std::string randomNumStr = std::to_string(randomNum);
    std::string timestampStr = std::to_string(timestamp);

    std::string result = timestampStr + randomNumStr;

    return result;
  }
}  // namespace micm