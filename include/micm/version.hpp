// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#pragma once

#ifdef __cplusplus
namespace micm {
extern "C" {
#endif

  inline const char* GetMicmVersion()
  {
    return "3.7.0";
  }
  inline unsigned GetMicmVersionMajor()
  {
    return 3;
  }
  inline unsigned GetMicmVersionMinor()
  {
    return 7+0;
  }
  inline unsigned GetMicmVersionPatch()
  {
    return 0+0;
  }
  inline unsigned GetMicmVersionTweak()
  {
    return +0;
  }

#ifdef __cplusplus
}  // extern "C"
}  // namespace micm
#endif
