// Copyright (C) 2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>

/// @brief Prefetches data to a specified locality level along the vector dimension for a vectorizable matrix
/// @param ptr Pointer to the first element to prefetch
/// @param size Total number of elements to prefetch
/// @param rw Read or write operation
/// @param locality Locality level (3 for L1, 2 for L2, 1 for L3)
#define PREFETCH_VECTOR(ptr, size, rw, locality) \
  { \
    const std::uint8_t* addr_byte = reinterpret_cast<const std::uint8_t*>(ptr); \
    for (std::size_t j = 0; j < size * sizeof(*ptr); j += prefetch::CACHE_LINE_SIZE) \
      { \
      __builtin_prefetch(addr_byte, rw, locality); \
      addr_byte += prefetch::CACHE_LINE_SIZE; \
    } \
  }

/// @brief Prefetches data to the L1 cache along the vector dimension for a vectorizable matrix
/// @param rw Read or write operation
/// @param ptr Pointer to the first element to prefetch
/// @param size Total number of elements to prefetch
#define PREFETCH_VECTOR_L1(rw, ptr, size) PREFETCH_VECTOR(ptr, size, rw, prefetch::L1_CACHE)

/// @brief Prefetches data to the L2 cache along the vector dimension for a vectorizable matrix
/// @param rw Read or write operation
/// @param ptr Pointer to the first element to prefetch
/// @param size Total number of elements to prefetch
#define PREFETCH_VECTOR_L2(rw, ptr, size) PREFETCH_VECTOR(ptr, size, rw, prefetch::L2_CACHE)

/// @brief Prefetches data to the L3 cache along the vector dimension for a vectorizable matrix
/// @param rw Read or write operation
/// @param ptr Pointer to the first element to prefetch
/// @param size Total number of elements to prefetch
#define PREFETCH_VECTOR_L3(rw, ptr, size) PREFETCH_VECTOR(ptr, size, rw, prefetch::L3_CACHE)

namespace micm
{
  namespace prefetch
  {
    /// @brief KiB in bytes
    constexpr std::size_t KiB = 1024;
    /// @brief MiB in bytes
    constexpr std::size_t MiB = 1024 * KiB;

    /// @brief L1 Cache locality hint
    constexpr int L1_CACHE = 3;
    /// @brief L2 Cache locality hint
    constexpr int L2_CACHE = 2;
    /// @brief L3 Cache locality hint
    constexpr int L3_CACHE = 1;

    /// @brief Prefetch read operation
    constexpr int READ = 0;
    /// @brief Prefetch write operation
    constexpr int WRITE = 1;
    
    // TODO The values below are machine-dependent.
    // We should add CMake options to set these values based on the target architecture.

    /// @brief Cache line size in bytes
    constexpr std::size_t CACHE_LINE_SIZE = 64;  // Assuming a cache line size of 64 bytes
    /// @brief L1 Cache size per core in bytes
    constexpr std::size_t L1_CACHE_SIZE = 32 * KiB;  // Assuming L1 cache size of 32 KiB
    /// @brief L2 Cache size per core in bytes
    constexpr std::size_t L2_CACHE_SIZE = 500 * KiB;  // Assuming L2 cache size of 500 KiB
    /// @brief L3 Cache size per core in bytes
    constexpr std::size_t L3_CACHE_SIZE = 8 * MiB;  // Assuming L3 cache size of 8 MiB

    

  }
}