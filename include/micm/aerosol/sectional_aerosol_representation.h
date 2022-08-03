/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_SECTIONAL_AEROSOL_REPRESENTATION_H
#define MICM_SECTIONAL_AEROSOL_REPRESENTATION_H

#include <micm/aerosol/aerosol_representation.h>
#include <micm/aerosol/section.h>

class SectionalAersolRepresentation: public AerosolRepresentation {
  private:
    const std::vector<Section> sections_;

  public:
};

#endif