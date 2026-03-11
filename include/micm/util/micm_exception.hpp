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

  /// @brief Base exception for all MICM errors. Catch this to handle any MICM error.
  struct MicmException : public std::runtime_error
  {
    MicmSeverity severity_;

    MicmException(MicmSeverity severity, const std::string& message)
        : std::runtime_error(message),
          severity_(severity)
    {
    }
  };

  /// @brief Typed MICM exception carrying a domain-specific error code.
  /// Catch this when you need to inspect or match the specific error code.
  template<typename ErrorCode>
  struct MicmCodedError : public MicmException
  {
    ErrorCode code_;

    MicmCodedError(ErrorCode code, MicmSeverity severity, const std::string& message)
        : MicmException(severity, message),
          code_(code)
    {
    }
  };
}  // namespace micm
