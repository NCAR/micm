// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/solver/temporary_variables.hpp>

namespace micm
{
  template<class DenseMatrixPolicy>
  class BackwardEulerDAETemporaryVariables : public TemporaryVariables
  {
   public:
    DenseMatrixPolicy Yn_;
    DenseMatrixPolicy forcing_;
    /// @brief Storage for algebraic constraint values
    std::vector<double> constraint_values_;

    BackwardEulerDAETemporaryVariables() = default;
    BackwardEulerDAETemporaryVariables(const BackwardEulerDAETemporaryVariables& other) = default;
    BackwardEulerDAETemporaryVariables(BackwardEulerDAETemporaryVariables&& other) = default;
    BackwardEulerDAETemporaryVariables& operator=(const BackwardEulerDAETemporaryVariables& other) = default;
    BackwardEulerDAETemporaryVariables& operator=(BackwardEulerDAETemporaryVariables&& other) = default;
    ~BackwardEulerDAETemporaryVariables() = default;

    BackwardEulerDAETemporaryVariables(const auto& state_parameters, const std::size_t number_of_grid_cells, const std::size_t number_of_constraints = 0)
        : Yn_(number_of_grid_cells, state_parameters.number_of_species_),
          forcing_(number_of_grid_cells, state_parameters.number_of_species_),
          constraint_values_(number_of_constraints, 0.0)
    {
    }
  };
}  // namespace micm