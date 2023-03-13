/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 */

#pragma once

#include <filesystem>
#include <iostream>
#include <micm/system/system.hpp>

namespace micm
{

  class JsonReaderPolicy
  {
   public:
    static std::unique_ptr<micm::System> ReadAndParse(std::filesystem::path path)
    {
      if (!std::filesystem::exists(path)) {
        std::string err_msg ="Configuration file at path " + path.string() + " does not exist\n";
        throw std::runtime_error(err_msg);
      }
      return std::make_unique<micm::System>(micm::System());
    }
  };

  template<class ConfigTypePolicy = JsonReaderPolicy>
  class SystemBuilder
  {
   public:
    std::unique_ptr<micm::System> Build(std::filesystem::path path)
    {
      return ConfigTypePolicy::ReadAndParse(path);
    }
  };

}  // namespace micm