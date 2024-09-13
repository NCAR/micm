// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/solver/temporary_variables.hpp>

namespace micm 
{
  template<class DenseMatrixPolicy>
  class BackwardEulerTemporaryVariables : public TemporaryVariables
  {
  public:
      DenseMatrixPolicy Yn_;
      DenseMatrixPolicy forcing_;

      BackwardEulerTemporaryVariables() = default;
      BackwardEulerTemporaryVariables(const BackwardEulerTemporaryVariables& other) = default;
      BackwardEulerTemporaryVariables(BackwardEulerTemporaryVariables&& other) = default;
      BackwardEulerTemporaryVariables& operator=(const BackwardEulerTemporaryVariables& other) = default;
      BackwardEulerTemporaryVariables& operator=(BackwardEulerTemporaryVariables&& other) = default;
      ~BackwardEulerTemporaryVariables() = default;

      BackwardEulerTemporaryVariables(const auto& state_parameters):
          Yn_(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_),
          forcing_(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_)
      {
      }
  };
} // namespace name