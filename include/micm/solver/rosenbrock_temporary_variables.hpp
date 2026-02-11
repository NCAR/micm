// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/solver/temporary_variables.hpp>

#include <vector>

namespace micm
{
  template<class DenseMatrixPolicy>
  class RosenbrockTemporaryVariables : public TemporaryVariables
  {
   public:
    DenseMatrixPolicy Ynew_;
    DenseMatrixPolicy initial_forcing_;
    std::vector<DenseMatrixPolicy> K_;
    DenseMatrixPolicy Yerror_;

    RosenbrockTemporaryVariables() = default;
    RosenbrockTemporaryVariables(const RosenbrockTemporaryVariables& other) = default;
    RosenbrockTemporaryVariables(RosenbrockTemporaryVariables&& other) = default;
    RosenbrockTemporaryVariables& operator=(const RosenbrockTemporaryVariables& other) = default;
    RosenbrockTemporaryVariables& operator=(RosenbrockTemporaryVariables&& other) = default;
    ~RosenbrockTemporaryVariables() = default;

    RosenbrockTemporaryVariables(
        const auto& state_parameters,
        const auto& solver_parameters,
        const std::size_t number_of_grid_cells)
        : Ynew_(
              number_of_grid_cells,
              state_parameters.constraints_replace_state_rows_
                  ? state_parameters.number_of_species_
                  : (state_parameters.number_of_species_ + state_parameters.number_of_constraints_)),
          initial_forcing_(
              number_of_grid_cells,
              state_parameters.constraints_replace_state_rows_
                  ? state_parameters.number_of_species_
                  : (state_parameters.number_of_species_ + state_parameters.number_of_constraints_)),
          Yerror_(
              number_of_grid_cells,
              state_parameters.constraints_replace_state_rows_
                  ? state_parameters.number_of_species_
                  : (state_parameters.number_of_species_ + state_parameters.number_of_constraints_))
    {
      const std::size_t matrix_columns = state_parameters.constraints_replace_state_rows_
                                             ? state_parameters.number_of_species_
                                             : (state_parameters.number_of_species_ + state_parameters.number_of_constraints_);
      K_.reserve(solver_parameters.stages_);
      for (std::size_t i = 0; i < solver_parameters.stages_; ++i)
        K_.emplace_back(number_of_grid_cells, matrix_columns);
    }
  };
}  // namespace micm
