/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 * 
 */

#pragma once

#include <micm/system/system.hpp>

#include <filesystem>

namespace micm {

class JsonReaderPolicy {

  public:
    static std::unique_ptr<micm::System> ReadAndParse(std::filesystem::path path) {
      return std::make_unique<micm::System>(micm::System());
    }

};

template<class ConfigTypePolicy = JsonReaderPolicy>
class SystemBuilder {
  public:

    std::unique_ptr<micm::System> Build(std::filesystem::path path) {
      return ConfigTypePolicy::ReadAndParse(path);
    }

};

}