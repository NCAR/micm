/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_MODAL_AEROSOL_REPRESENTATION_H
#define MICM_MODAL_AEROSOL_REPRESENTATION_H

#include <micm/system/aerosol_representation.h>
#include <micm/system/mode.h>

class ModalAerosolRepresentation: public AerosolRepresentation {
  private:
    const std::vector<Mode> modes_;

  public:
    virtual void GenerateJITForSpecies(JIT jit, const std::function<void(JIT, std::map<std::string, Species>)>& f) = 0;
};

#endif