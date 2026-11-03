// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>

namespace micm
{
  /// @brief Generates a version 4 random Universally Unique Identifier
  ///        Used to build parameter names that are guaranteed unique across models
  /// @return UUID string in the form xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx, where y is one of 8, 9, a, b
  inline std::string GenerateUuid()
  {
    // One engine per thread, seeded once from several random_device draws. Re-seeding on every
    // call from a single 32-bit draw caps the entropy at 32 bits (less where random_device is weak)
    // and produced duplicate UUIDs in practice.
    thread_local std::mt19937_64 gen = []
    {
      std::random_device rd;
      std::seed_seq seq{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
      return std::mt19937_64(seq);
    }();
    std::uniform_int_distribution<std::uint64_t> dis;

    std::uint64_t high = dis(gen);
    std::uint64_t low = dis(gen);
    high = (high & 0xFFFFFFFFFFFF0FFFULL) | 0x0000000000004000ULL;  // UUID version 4
    low = (low & 0x3FFFFFFFFFFFFFFFULL) | 0x8000000000000000ULL;    // Variant bits (10xx)

    std::ostringstream ss;
    ss << std::hex << std::setfill('0') << std::setw(8) << (high >> 32) << "-" << std::setw(4) << ((high >> 16) & 0xFFFF)
       << "-" << std::setw(4) << (high & 0xFFFF) << "-" << std::setw(4) << (low >> 48) << "-" << std::setw(12)
       << (low & 0xFFFFFFFFFFFFULL);
    return ss.str();
  }
}  // namespace micm
