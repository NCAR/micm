// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdint>

#if defined(MICM_USE_SINGLE) && defined(MICM_USE_DOUBLE)
  #error "MICM_USE_SINGLE and MICM_USE_DOUBLE are mutually exclusive"
#endif

#ifdef MICM_USE_KOKKOS
  #include <Kokkos_Macros.hpp>
  #define MICM_LAMBDA                 KOKKOS_LAMBDA
  #define MICM_DEVICE_FUNCTION        KOKKOS_INLINE_FUNCTION
  #define MICM_INLINE_DEVICE_FUNCTION KOKKOS_INLINE_FUNCTION
#else
  #define MICM_LAMBDA [=]
  #define MICM_DEVICE_FUNCTION
  #define MICM_INLINE_DEVICE_FUNCTION inline
#endif

#ifndef MICM_DEFAULT_VECTOR_SIZE
  #define MICM_DEFAULT_VECTOR_SIZE 4
#endif

namespace micm
{

// Real-precision typedef.  Single precision is opt-in, so that this header alone -- included
// without the micm cmake target to supply any definitions -- yields the same type a default cmake
// build does.  Both macros are propagated by the target: cmake -DMICM_USE_DOUBLE=OFF (the
// non-default) defines MICM_USE_SINGLE, and MICM_USE_DOUBLE=ON (the default) defines
// MICM_USE_DOUBLE, so downstream code can branch on either without inferring precision from the
// absence of a macro.
#ifdef MICM_USE_SINGLE
  using Real = float;
#else
  using Real = double;
#endif

  using Index = std::size_t;

  // Boolean stored as a single byte.  Used in place of Bool wherever a contiguous array of flags is
  // needed: std::vector<Bool> is a bit-packed specialization, so it offers no data() pointer to hand
  // to a CUDA memcpy or to index from device code.
  using Bool = std::uint8_t;

  /// @brief A device-compatible struct for holding three indices
  struct IndexTrio
  {
    Index first_;
    Index second_;
    Index third_;
  };

  /// @brief A device-compatible struct for holding pairs of indices
  struct IndexPair
  {
    Index first_;
    Index second_;
  };

}  // namespace micm
