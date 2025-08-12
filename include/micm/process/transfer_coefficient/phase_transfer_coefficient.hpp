// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/conditions.hpp>

#include <memory>

namespace micm
{
  // TODO - This is based on the implementation in cam-acom-dev
  //
  /// @brief van 't Hoff equation parameters
  /// Defines the temperature dependence of Henry's Law constants using:
  ///     H(T) = A * exp(-B * (1/T - 1/T_ref))
  /// where:
  ///     - A and B are the parameters,
  ///     - T is the current temperature,
  ///     - T_ref is the reference temperature.
  struct Van_t_HoffParams
  {
    double C = 0.0;  // mol L-1 || mol2 L-2 || mol L-1 atm-1
    double T = 0.0;  // K

    Van_t_HoffParams() = default;

    Van_t_HoffParams(double C, double T)
        : C(C),
          T(T)
    {
    }
  };

  // TODO - class PhaseTransferCoefficient is a placeholder and will be further developed.
  //        issue: https://github.com/NCAR/micm/issues/811
  class PhaseTransferCoefficient : public TransferCoefficient
  {
   public:
    PhaseTransferCoefficient() = default;

    virtual ~PhaseTransferCoefficient() = default;

    virtual std::unique_ptr<TransferCoefficient> Clone() const override
    {
      return std::make_unique<PhaseTransferCoefficient>(*this);
    }

    virtual double Calculate() const override
    {
      return 1.0;
    }

    virtual double Calculate(const Conditions& conditions) const override
    {
      return 1.0;
    }
  };
}  // namespace micm