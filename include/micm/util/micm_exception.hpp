// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>

#include <stdexcept>
#include <string>
#include <string_view>

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
    std::string_view category_;
    int code_;

    MicmException(MicmSeverity severity, std::string_view category, int code, const std::string& message)
        : std::runtime_error(message),
          severity_(severity),
          category_(category),
          code_(code)
    {
    }
  };

}  // namespace micm
