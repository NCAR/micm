/* Copyright (C) 2023-2024 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/conditions.hpp>
#include <micm/util/error.hpp>

#include <cstddef>
#include <iterator>
#include <memory>
#include <system_error>
#include <vector>

enum class MicmRateConstantErrc
{
  MissingArgumentsForSurfaceRateConstant = MICM_RATE_CONSTANT_ERROR_CODE_MISSING_ARGUMENTS_FOR_SURFACE_RATE_CONSTANT,
  MissingArgumentsForUserDefinedRateConstant = MICM_RATE_CONSTANT_ERROR_CODE_MISSING_ARGUMENTS_FOR_USER_DEFINED_RATE_CONSTANT
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmRateConstantErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class RateConstantErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_RATE_CONSTANT;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmRateConstantErrc>(ev))
      {
        case MicmRateConstantErrc::MissingArgumentsForSurfaceRateConstant: return "Missing arguments for surface reaction";
        case MicmRateConstantErrc::MissingArgumentsForUserDefinedRateConstant:
          return "Missing arguments for user-defined rate constant";
        default: return "Unknown error";
      }
    }
  };

  const RateConstantErrorCategory rateConstantErrorCategory{};
}  // namespace

std::error_code make_error_code(MicmRateConstantErrc e)
{
  return { static_cast<int>(e), rateConstantErrorCategory };
}

namespace micm
{
  /// @brief A base class for any type of rate constant
  class RateConstant
  {
   public:
    /// @brief Virtual destructor
    virtual ~RateConstant(){};

    /// @brief Deep copy
    virtual std::unique_ptr<RateConstant> clone() const = 0;

    /// @brief Returns a set of labels for user-defined rate constant parameters
    /// @return Vector of custom parameter labels
    virtual std::vector<std::string> CustomParameters() const
    {
      return std::vector<std::string>{};
    }

    /// @brief Returns the number of custom parameters
    /// @return Number of custom parameters
    virtual std::size_t SizeCustomParameters() const
    {
      return 0;
    }

    /// @brief Calculate the rate constant for a set of conditions
    /// @param conditions The current environmental conditions of the chemical system
    /// @return The reaction rate constant
    virtual double calculate(const Conditions& conditions) const
    {
      return 0;
    }

    /// @brief Calculate the rate constant for a set of conditions
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User defined rate constant parameters
    /// @return The reaction rate constant
    virtual double calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const
    {
      return 0;
    }
  };

}  // namespace micm