/* Copyright (C) 2023-2024 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
// clang-format off
#pragma once

#ifdef __cplusplus
namespace micm {
extern "C" {
#endif

  const char* getMicmVersion()
  {
    return "3.5.0";
  }
  unsigned getMicmVersionMajor()
  {
    return 3;
  }
  unsigned getMicmVersionMinor()
  {
    return 5+0;
  }
  unsigned getMicmVersionPatch()
  {
    return 0+0;
  }
  unsigned getMicmVersionTweak()
  {
    return +0;
  }

#ifdef __cplusplus
}  // extern "C"
}  // namespace micm
#endif
