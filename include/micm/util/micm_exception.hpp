// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>

#include <optional>
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

  struct MicmException : public std::runtime_error
  {
    MicmSeverity severity_;
    const char* category_;
    int code_;

    MicmException(MicmSeverity severity, const char* category, int code, const std::string& message)
        : std::runtime_error(message),
          severity_(severity),
          category_(category),
          code_(code)
    {
    }
  };

  /// @brief Wraps a successfully constructed value alongside an optional warning.
  ///        Use with TryBuild() to receive the built object even when non-fatal
  ///        issues (e.g. unused species) are detected.
  template<typename T>
  struct MicmResult
  {
    T value_;
    std::optional<MicmException> warning_;
  };

}  // namespace micm
