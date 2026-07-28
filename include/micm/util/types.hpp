// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdint>

namespace micm
{

// Real-precision typedef.  Default ON via the MICM_USE_DOUBLE cmake option.
#ifdef MICM_USE_DOUBLE
  using Real = double;
#else
  using Real = float;
#endif

  using Index = std::size_t;

  // Boolean stored as a single byte.  Used in place of bool wherever a contiguous array of flags is
  // needed: std::vector<bool> is a bit-packed specialization, so it offers no data() pointer to hand
  // to a CUDA memcpy or to index from device code.
  using Bool = std::uint8_t;

}  // namespace micm
