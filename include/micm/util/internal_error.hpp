// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>
#include <micm/util/micm_exception.hpp>

#include <string>

#define INTERNAL_ERROR(msg) micm::ThrowInternalError(MicmInternalErrc::General, __FILE__, __LINE__, msg);

enum class MicmInternalErrc
{
  General = MICM_INTERNAL_ERROR_CODE_GENERAL,
  Cuda = MICM_INTERNAL_ERROR_CODE_CUDA,
  Cublas = MICM_INTERNAL_ERROR_CODE_CUBLAS
};

namespace micm
{
  inline void ThrowInternalError(MicmInternalErrc e, const char* file, int line, const char* msg)
  {
    std::string message = std::string("Please file a bug report at https://github.com/NCAR/micm. Error detail: (") +
                          file + ":" + std::to_string(line) + ") " + msg;
    throw MicmException<MicmInternalErrc>(e, MicmSeverity::Critical, message);
  }
}  // namespace micm
