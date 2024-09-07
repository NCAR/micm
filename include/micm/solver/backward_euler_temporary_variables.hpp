// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace micm 
{
  template<class DenseMatrixPolicy>
  struct BackwardEulerTemporaryVariables
  {
  public:
      BackwardEulerTemporaryVariables() = default;
      BackwardEulerTemporaryVariables(const BackwardEulerTemporaryVariables& other) = delete;
      BackwardEulerTemporaryVariables(BackwardEulerTemporaryVariables&& other) = default;
      BackwardEulerTemporaryVariables& operator=(const BackwardEulerTemporaryVariables& other) = delete;
      BackwardEulerTemporaryVariables& operator=(BackwardEulerTemporaryVariables&& other) = default;
      ~BackwardEulerTemporaryVariables() = default;

      BackwardEulerTemporaryVariables(const auto& state, const auto& parameters)
      { }
  };
} // namespace name