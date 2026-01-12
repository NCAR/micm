// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
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

    BackwardEulerTemporaryVariables(const auto& state_parameters, const std::size_t number_of_grid_cells)
        : Yn_(number_of_grid_cells, state_parameters.number_of_species_),
          forcing_(number_of_grid_cells, state_parameters.number_of_species_)
    {
    }
  };
}  // namespace micm