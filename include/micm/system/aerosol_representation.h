/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_AEROSOL_REPRESENTATION_H
#define MICM_AEROSOL_REPRESENTATION_H

#include <micm/jit.h>
#include <micm/system/species.h>

#include <cstddef>
#include <map>
#include <string>

class AerosolRepresentation {
  private:

  public:
    std::size_t Size();
    virtual void GenerateJITForSpecies(JIT jit, const std::function<void(JIT, std::map<std::string, Species>)>& f) = 0;
};

#endif