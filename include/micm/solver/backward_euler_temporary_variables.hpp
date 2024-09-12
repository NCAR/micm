// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace micm 
{
  template<class DenseMatrixPolicy>
  class BackwardEulerTemporaryVariables
  {
  public:
      DenseMatrixPolicy Yn_;
      DenseMatrixPolicy forcing_;

      BackwardEulerTemporaryVariables() = default;
      BackwardEulerTemporaryVariables(const BackwardEulerTemporaryVariables& other) = delete;
      BackwardEulerTemporaryVariables(BackwardEulerTemporaryVariables&& other) = default;
      BackwardEulerTemporaryVariables& operator=(const BackwardEulerTemporaryVariables& other) = delete;
      BackwardEulerTemporaryVariables& operator=(BackwardEulerTemporaryVariables&& other) = default;
      ~BackwardEulerTemporaryVariables() = default;

      BackwardEulerTemporaryVariables(const auto& state_parameters):
          Yn_(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_),
          forcing_(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_)
      {
      }
  };
} // namespace name