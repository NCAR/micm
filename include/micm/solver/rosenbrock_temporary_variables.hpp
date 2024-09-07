// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

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

        RosenbrockTemporaryVariables() = default;
        RosenbrockTemporaryVariables(const RosenbrockTemporaryVariables& other) = delete;
        RosenbrockTemporaryVariables(RosenbrockTemporaryVariables&& other) = default;
        RosenbrockTemporaryVariables& operator=(const RosenbrockTemporaryVariables& other) = delete;
        RosenbrockTemporaryVariables& operator=(RosenbrockTemporaryVariables&& other) = default;
        ~RosenbrockTemporaryVariables() = default;
        
        RosenbrockTemporaryVariables(const auto& state, const auto& parameters) :
            Ynew_(state.variables_.NumRows(), state.variables_.NumColumns()),
            initial_forcing_(state.variables_.NumRows(), state.variables_.NumColumns()),
            Yerror_(state.variables_.NumRows(), state.variables_.NumColumns())
        {
            K_.reserve(parameters.stages_);
            for (std::size_t i = 0; i < parameters.stages_; ++i)
                K_.emplace_back(state.variables_.NumRows(), state.variables_.NumColumns());
        }
    };
}  // namespace micm