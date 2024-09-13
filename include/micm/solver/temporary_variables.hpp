// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace micm
{
    /// @brief This is the base class for temporary variables; currently it is empty and will be expanded by a specific solver later.
    class TemporaryVariables
    {
    public:
        TemporaryVariables() = default;
        TemporaryVariables(const TemporaryVariables& other) = delete;
        TemporaryVariables(TemporaryVariables&& other) = default;
        TemporaryVariables& operator=(const TemporaryVariables& other) = delete;
        TemporaryVariables& operator=(TemporaryVariables&& other) = default;
        ~TemporaryVariables() = default;
    };
}  // namespace micm