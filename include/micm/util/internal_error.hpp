// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>
#include <micm/util/micm_exception.hpp>

#include <string>

#define INTERNAL_ERROR(msg) micm::ThrowInternalError(__FILE__, __LINE__, msg);

namespace micm
{
  inline void ThrowInternalError(const char* file, int line, const char* msg)
  {
    std::string message = std::string("Please file a bug report at https://github.com/NCAR/micm. Error detail: (") +
                          file + ":" + std::to_string(line) + ") " + msg;
    throw MicmException(MicmSeverity::Critical, MICM_ERROR_CATEGORY_INTERNAL, MICM_INTERNAL_ERROR_CODE_GENERAL, message);
  }
}  // namespace micm
