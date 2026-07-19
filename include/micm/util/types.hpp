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

}  // namespace micm
