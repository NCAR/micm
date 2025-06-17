// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant.hpp>

#include <string>
#include <functional>

namespace micm
{
  struct LambdaRateConstantParameters
  {
    /// @brief Label for the reaction used to identify user-defined parameters
    std::string label_;
    /// @brief Lambda function for calculating the rate constant
    std::function<double(const Conditions&)> lambda_function_;
  };

  /// @brief A lambda-backed rate constant
  class LambdaRateConstant : public RateConstant
  {
   public:
    LambdaRateConstantParameters parameters_;

    /// @brief Default constructor.
    LambdaRateConstant();

    /// @brief
    /// @param parameters The data needed to build this class
    LambdaRateConstant(const LambdaRateConstantParameters& parameters);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> Clone() const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double Calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return A rate constant based off of the conditions in the system
    double Calculate(const Conditions& conditions) const override;
  };

  inline LambdaRateConstant::LambdaRateConstant()
      : parameters_()
  {
  }

  inline LambdaRateConstant::LambdaRateConstant(const LambdaRateConstantParameters& parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<RateConstant> LambdaRateConstant::Clone() const
  {
    return std::unique_ptr<RateConstant>{ new LambdaRateConstant{ *this } };
  }

  inline double LambdaRateConstant::Calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    return parameters_.lambda_function_(conditions);
  }

  inline double LambdaRateConstant::Calculate(const Conditions& conditions) const
  {
    return parameters_.lambda_function_(conditions);
  }
}