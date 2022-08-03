/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
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
};

#endif