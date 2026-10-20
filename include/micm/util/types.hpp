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

// MICM_NO_UNROLL stops GCC from unrolling (and completely peeling) the loop that follows it.
// micm::VectorLoop uses it. The loop is still vectorized.
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_LLVM_COMPILER) && !defined(__NVCOMPILER) && \
    !defined(__CUDACC__)
  #define MICM_NO_UNROLL _Pragma("GCC unroll 1")
#else
  #define MICM_NO_UNROLL
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

  /// @brief Calls body(i) for each i in [0, L), for the L grid cells of a vector-ordered group
  ///
  /// GCC completely peels a loop of up to 16 vector iterations (max-completely-peel-times). For vector
  /// lengths 16 to 32 on aarch64, the peeled code ran 3 to 4 times slower than a loop, so for L >= 16 this
  /// keeps the loop. For L < 16, the complete unrolling is faster, so the loop is a plain loop.
  /// A pragma cannot depend on a template argument, so `if constexpr` selects the loop.
  template<Index L, typename Body>
  [[gnu::always_inline]] inline void VectorLoop(Body&& body)
  {
    if constexpr (L >= 16)
    {
      MICM_NO_UNROLL
      for (Index i = 0; i < L; ++i)
      {
        body(i);
      }
    }
    else
    {
      for (Index i = 0; i < L; ++i)
      {
        body(i);
      }
    }
  }

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

  /// @brief A 32-bit index for the large index lists of the solver kernels
  ///        The index lists of the in-place Mozart LU decomposition hold offsets inside one group of the sparse
  ///        matrix, which fit in 32 bits for any practical mechanism. With half the size, a kernel at a small vector
  ///        length loads half the index bytes for each multiply-add.
  using CompactIndex = std::uint32_t;

  /// @brief A device-compatible struct for holding three 32-bit indices
  struct CompactIndexTrio
  {
    CompactIndex first_;
    CompactIndex second_;
    CompactIndex third_;
  };

  /// @brief A device-compatible struct for holding pairs of 32-bit indices
  struct CompactIndexPair
  {
    CompactIndex first_;
    CompactIndex second_;
  };

}  // namespace micm
