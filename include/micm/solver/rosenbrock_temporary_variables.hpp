// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/solver/rosenbrock_solver_parameters.hpp>

namespace micm
{
    template<class DenseMatrixPolicy>
    struct RosenbrockTemporaryVariables
    {
    public:
        DenseMatrixPolicy Ynew_;
        DenseMatrixPolicy initial_forcing_;
        std::vector<DenseMatrixPolicy> K_;
        DenseMatrixPolicy Yerror_;

        RosenbrockTemporaryVariables() = delete;
        RosenbrockTemporaryVariables(const RosenbrockTemporaryVariables& other) = delete;
        RosenbrockTemporaryVariables(RosenbrockTemporaryVariables&& other) = default;
        RosenbrockTemporaryVariables& operator=(const RosenbrockTemporaryVariables& other) = delete;
        RosenbrockTemporaryVariables& operator=(RosenbrockTemporaryVariables&& other) = default;
        ~RosenbrockTemporaryVariables() = default;
        
        RosenbrockTemporaryVariables(const auto& state, const RosenbrockSolverParameters& parameters) :
            Ynew_(state.variables_.NumRows(), state.variables_.NumColumns()),
            initial_forcing_(state.variables_.NumRows(), state.variables_.NumColumns()),
            Yerror_(state.variables_.NumRows(), state.variables_.NumColumns())
        {
            K.reserve(parameters.stages_);
            for (std::size_t i = 0; i < parameters.stages_; ++i)
                K.emplace_back(num_rows, num_cols);
        }
    };
}  // namespace micm