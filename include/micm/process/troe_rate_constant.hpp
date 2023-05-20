/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cmath>
#include <micm/process/rate_constant.hpp>

namespace micm
{

  /**
   * @brief A Troe rate constant
   *
   */
  class TroeRateConstant : public RateConstant
  {
   private:
    /// @brief // TODO:
    const double k0_A_;
    /// @brief // TODO:
    const double k0_B_;
    /// @brief // TODO:
    const double k0_C_;
    /// @brief // TODO:
    const double kinf_A_;
    /// @brief // TODO:
    const double kinf_B_;
    /// @brief // TODO:
    const double kinf_C_;
    /// @brief // TODO:
    const double Fc_;
    /// @brief // TODO:
    const double N_;

   public:
    /// @brief Default constructor
    TroeRateConstant();

    /// @brief Calculate the rate constant
    /// @param state The current state of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double calculate(const State& state, std::vector<double>::const_iterator custom_parameters) const override;

    /// @brief Deep copy
    std::unique_ptr<RateConstant> clone() const override;

    /// @brief Calculate the rate constant
    /// @param temperature Temperature in [K]
    /// @param air_number_density Number density in [# cm-3]
    /// @return
    double calculate(const double& temperature, const double& air_number_density) const;
  };

  inline TroeRateConstant::TroeRateConstant()
      : k0_A_(),
        k0_B_(),
        k0_C_(),
        kinf_A_(),
        kinf_B_(),
        kinf_C_(),
        Fc_(),
        N_()
  {
  }

  inline std::unique_ptr<RateConstant> TroeRateConstant::clone() const {
    return std::unique_ptr<RateConstant>{new TroeRateConstant{*this}};
  }

  inline double TroeRateConstant::calculate(const State& state, std::vector<double>::const_iterator custom_parameters) const
  {
    return calculate(state.temperature_, state.air_density_);
  }

  inline double TroeRateConstant::calculate(const double& temperature, const double& air_number_density) const
  {
    double k0 = this->k0_A_ * std::exp(this->k0_C_ / temperature) * pow(temperature / 300.0, this->k0_B_);
    double kinf = this->kinf_A_ * std::exp(this->kinf_C_ / temperature) * pow(temperature / 300.0, this->kinf_B_);

    return k0 * air_number_density / (1.0 + k0 * air_number_density / kinf) *
           pow(this->Fc_, 1.0 / (1.0 + 1.0 / this->N_ * pow(log10(k0 * air_number_density / kinf), 2)));
  }

}  // namespace micm