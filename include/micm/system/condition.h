/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_CONDITION_H
#define MICM_CONDITION_H

#include <string>

class Condition {
  private:
    const std::string name_;
    const std::string units_;

  public:
    Condition(const std::string& name, const std::string& units);
};

#endif