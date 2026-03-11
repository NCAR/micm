// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stdexcept>
#include <string>

namespace micm
{
  enum class MicmSeverity
  {
    Warning,
    Error,
    Critical
  };

  template<typename ErrorCode>
  struct MicmException : public std::runtime_error
  {
    ErrorCode code_;
    MicmSeverity severity_;

    MicmException(ErrorCode code, MicmSeverity severity, const std::string& message)
        : std::runtime_error(message),
          code_(code),
          severity_(severity)
    {
    }
  };
}  // namespace micm
