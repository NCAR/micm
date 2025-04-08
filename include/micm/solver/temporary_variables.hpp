// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace micm
{
  /// @brief This is the base class for temporary variables; currently it is empty and will be expanded by a specific solver
  /// later.
  class TemporaryVariables
  {
   public:
    TemporaryVariables() = default;
    TemporaryVariables(const TemporaryVariables& other) = default;
    TemporaryVariables(TemporaryVariables&& other) = default;
    TemporaryVariables& operator=(const TemporaryVariables& other) = default;
    TemporaryVariables& operator=(TemporaryVariables&& other) = default;
    virtual ~TemporaryVariables() = default;
  };
}  // namespace micm