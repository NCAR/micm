/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_SECTION_H
#define MICM_SECTION_H

#include <micm/system/phase.h>

class Section {
  private:
    const std::vector<const Phase&> phases_;

  public:
    std::size_t Size();
};

#endif