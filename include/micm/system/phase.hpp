/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/species.hpp>
#include <vector>

namespace micm
{

  class Phase
  {
   private:
    const std::vector<Species> species_;

   public:
    Phase();
    Phase(const Phase& other);
    Phase(std::vector<Species> species);

    Phase& operator=(const Phase& other);

    size_t Size();
  };

}  // namespace micm