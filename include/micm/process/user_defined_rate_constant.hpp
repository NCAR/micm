// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant.hpp>
#include <string>

namespace micm
{
  struct UserDefinedRateConstantParameters
  {
    /// @brief Label for the reaction used to identify user-defined parameters
    std::string label_;
    /// @brief Scaling factor to apply to user-provided rate constants
    double scaling_factor_{ 1.0 };
  };

  /// @brief A photolysis rate constant
  class UserDefinedRateConstant : public RateConstant
  {
   public:
    UserDefinedRateConstantParameters parameters_;

    /// @brief Default constructor.
    UserDefinedRateConstant();

    /// @brief
    /// @param parameters The data needed to build this class
    UserDefinedRateConstant(const UserDefinedRateConstantParameters& parameters);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> clone() const override;

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
    double calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return A rate constant based off of the conditions in the system
    double calculate(const Conditions& conditions) const override;
  };

  inline UserDefinedRateConstant::UserDefinedRateConstant()
      : parameters_()
  {
  }

  inline UserDefinedRateConstant::UserDefinedRateConstant(const UserDefinedRateConstantParameters& parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<RateConstant> UserDefinedRateConstant::clone() const
  {
    return std::unique_ptr<RateConstant>{ new UserDefinedRateConstant{ *this } };
  }

  inline double UserDefinedRateConstant::calculate(const Conditions& conditions) const
  {
    throw std::runtime_error(
        "User defined rate constants must be supplied with custom rate parameters using the alternative calculate function");
  }

  inline double UserDefinedRateConstant::calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
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