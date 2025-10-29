// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant/rate_constant.hpp>

#include <string>
#include <functional>

namespace micm
{
  struct UserDefinedRateConstantParameters
  {
    /// @brief Label for the reaction used to identify user-defined parameters
    std::string label_;
    /// @brief Scaling factor to apply to user-provided rate constants
    double scaling_factor_{ 1.0 };
    /// @brief Optional function to compute the rate constant
    /// Takes temperature, pressure, and air density as parameters
    /// and returns the computed rate constant.
    /// If not provided, a user-supplied constant value is used instead.
    std::function<double(double, double, double)> parameterize_{ nullptr };
  };

  /// @brief A photolysis rate constant
  class UserDefinedRateConstant : public RateConstant
  {
    /// @brief A function that if provided will be used to return a rate constant
    std::function<double(double, double, double)> parameterize_{ nullptr };

   public:
    UserDefinedRateConstantParameters parameters_;

    /// @brief Default constructor.
    UserDefinedRateConstant();

    /// @brief
    /// @param parameters The data needed to build this class
    UserDefinedRateConstant(const UserDefinedRateConstantParameters& parameters);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> Clone() const override;

    /// @brief Returns a label for the user-defined rate constant parameter
    /// @return Rate constant label
    std::vector<std::string> CustomParameters() const override;

    /// @brief Returns the number of custom parameters
    /// @return Number of custom parameters
    std::size_t SizeCustomParameters() const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double Calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;

    bool IsParameterized() const override
    {
      return parameterize_ != nullptr;
    }
  };

  inline UserDefinedRateConstant::UserDefinedRateConstant()
      : parameters_()
      , parameterize_(nullptr)
  {
  }

  inline UserDefinedRateConstant::UserDefinedRateConstant(const UserDefinedRateConstantParameters& parameters)
      : parameters_(parameters)
      , parameterize_(parameters.parameterize_)
  {
  }

  inline std::unique_ptr<RateConstant> UserDefinedRateConstant::Clone() const
  {
    return std::make_unique<UserDefinedRateConstant>(*this);
  }

  inline double UserDefinedRateConstant::Calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    if (parameterize_)
    {
      return parameterize_(conditions.temperature_, conditions.pressure_, conditions.air_density_) * parameters_.scaling_factor_;
    }
    return (double)*custom_parameters * parameters_.scaling_factor_;
  }

  inline std::vector<std::string> UserDefinedRateConstant::CustomParameters() const
  {
    return std::vector<std::string>{ parameters_.label_ };
  }

  inline std::size_t UserDefinedRateConstant::SizeCustomParameters() const
  {
    return 1;
  }
}  // namespace micm