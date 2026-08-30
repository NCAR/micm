// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/solver/temporary_variables.hpp>
#include <micm/util/types.hpp>

#include <vector>

namespace micm
{
  template<class DenseMatrixPolicy>
  class RosenbrockTemporaryVariables : public TemporaryVariables
  {
    template<class U>
    using Scalar = typename DenseMatrixPolicy::template ScalarType<U>;

   public:
    DenseMatrixPolicy Ynew_;
    DenseMatrixPolicy initial_forcing_;
    std::vector<DenseMatrixPolicy> K_;
    DenseMatrixPolicy Yerror_;
    Scalar<Real> current_c_over_h_;
    Scalar<Real> error_;
    Scalar<Real> max_residual_;
    Scalar<Real> max_correction_;
    Scalar<Bool> nan_detected_;
    Scalar<Bool> inf_detected_;
    // Line-search step length for constraint initialization. This has to be a device-resident
    // scalar handle rather than a plain Real: MICM_LAMBDA captures by value, so a Real named
    // inside a Function built once outside the backtrack loop would be frozen at construction.
    Scalar<Real> constraint_init_step_;

    RosenbrockTemporaryVariables() = default;
    RosenbrockTemporaryVariables(const RosenbrockTemporaryVariables& other) = default;
    RosenbrockTemporaryVariables(RosenbrockTemporaryVariables&& other) = default;
    RosenbrockTemporaryVariables& operator=(const RosenbrockTemporaryVariables& other) = default;
    RosenbrockTemporaryVariables& operator=(RosenbrockTemporaryVariables&& other) = default;
    ~RosenbrockTemporaryVariables() override = default;

    std::unique_ptr<TemporaryVariables> Clone() const override
    {
      return std::make_unique<RosenbrockTemporaryVariables>(*this);
    }

    RosenbrockTemporaryVariables(
        const auto& state_parameters,
        const auto& solver_parameters,
        const Index number_of_grid_cells)
        : Ynew_(number_of_grid_cells, state_parameters.number_of_species_),
          initial_forcing_(number_of_grid_cells, state_parameters.number_of_species_),
          Yerror_(number_of_grid_cells, state_parameters.number_of_species_)
    {
      K_.reserve(solver_parameters.stages_);
      for (Index i = 0; i < solver_parameters.stages_; ++i)
      {
        K_.emplace_back(number_of_grid_cells, state_parameters.number_of_species_);
      }
    }
  };
}  // namespace micm
