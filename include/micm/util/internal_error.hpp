/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/util/error.hpp>

#include <system_error>

#define INTERNAL_ERROR(msg) micm::ThrowInternalError(MicmInternalErrc::General, __FILE__, __LINE__, msg);

enum class MicmInternalErrc
{
  General = MICM_INTERNAL_ERROR_CODE_GENERAL,
  Cuda = MICM_INTERNAL_ERROR_CODE_CUDA,
  Cublas = MICM_INTERNAL_ERROR_CODE_CUBLAS
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmInternalErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class MicmInternalErrorCategory : public std::error_category
  {
   public:
    const char *name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_INTERNAL;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmInternalErrc>(ev))
      {
        case MicmInternalErrc::General: return "Internal error";
        case MicmInternalErrc::Cuda: return "CUDA error";
        case MicmInternalErrc::Cublas: return "cuBLAS error";
        default: return "Unknown error";
      }
    }
  };

  const MicmInternalErrorCategory micmInternalErrorCategory{};
}  // namespace

inline std::error_code make_error_code(MicmInternalErrc e)
{
  return { static_cast<int>(e), micmInternalErrorCategory };
}

namespace micm
{
  inline void ThrowInternalError(MicmInternalErrc e, const char *file, int line, const char *msg)
  {
    std::string message = std::string("Please file a bug report at https://github.com/NCAR/micm. Error detail: (") + file +
                          ":" + std::to_string(line) + ") " + msg;
    throw std::system_error(make_error_code(e), message);
  }
}  // namespace micm