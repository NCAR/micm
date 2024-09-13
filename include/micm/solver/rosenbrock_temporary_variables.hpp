// Copyright (C) 2023-2024 National Center for Atmospheric Research
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

        RosenbrockTemporaryVariables(const auto& state_parameters, const auto& solver_parameters):
            Ynew_(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_),
            initial_forcing_(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_),
            Yerror_(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_)
        {
            K_.reserve(solver_parameters.stages_);
            for (std::size_t i = 0; i < solver_parameters.stages_; ++i)
                K_.emplace_back(state_parameters.number_of_grid_cells_, state_parameters.number_of_species_);
        }
    };
}  // namespace micm