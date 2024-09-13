// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace micm
{
    class EmptyTemporaryVariables
    {
    public:
        EmptyTemporaryVariables() = default;
        EmptyTemporaryVariables(const EmptyTemporaryVariables& other) = delete;
        EmptyTemporaryVariables(EmptyTemporaryVariables&& other) = default;
        EmptyTemporaryVariables& operator=(const EmptyTemporaryVariables& other) = delete;
        EmptyTemporaryVariables& operator=(EmptyTemporaryVariables&& other) = default;
        ~EmptyTemporaryVariables() = default;
    };
}  // namespace micm