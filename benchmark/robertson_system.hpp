// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Robertson stiff problem definition, shared by the full-ODE and DAE benchmark
// solver builds so both run provably identical chemistry.
//
//   r1: A      -> B         rate = k1 * [A]
//   r2: 2B     -> B + C     rate = k2 * [B]^2
//   r3: B + C  -> A + C     rate = k3 * [B] * [C]
#pragma once

#include <micm/CPU.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace robertson
{
  /// Default Robertson rate constants (Hairer & Wanner II, p.3).
  inline constexpr double K1_DEFAULT = 0.04;
  inline constexpr double K2_DEFAULT = 3.0e7;  // the stiffness knob
  inline constexpr double K3_DEFAULT = 1.0e4;

  /// The species and reactions of the Robertson system. The gas phase must be
  /// kept alive for the lifetime of the solver, so it is returned alongside the
  /// processes.
  struct System
  {
    micm::Phase gas_phase;
    std::vector<micm::Process> processes;
  };

  /// Build the Robertson reaction system. Rate constants are supplied at solve
  /// time via SetCustomRateParameter("r1"/"r2"/"r3", ...), not baked in here.
  inline System MakeSystem()
  {
    auto a = micm::Species("A");
    auto b = micm::Species("B");
    auto c = micm::Species("C");

    micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

    micm::Process r1 = micm::ChemicalReactionBuilder()
                           .SetReactants({ a })
                           .SetProducts({ micm::StoichSpecies(b, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r1" })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process r2 = micm::ChemicalReactionBuilder()
                           .SetReactants({ b, b })
                           .SetProducts({ micm::StoichSpecies(b, 1), micm::StoichSpecies(c, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r2" })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process r3 = micm::ChemicalReactionBuilder()
                           .SetReactants({ b, c })
                           .SetProducts({ micm::StoichSpecies(a, 1), micm::StoichSpecies(c, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r3" })
                           .SetPhase(gas_phase)
                           .Build();

    return System{ gas_phase, { r1, r2, r3 } };
  }

  /// Consistent quasi-steady-state value of B given A and C (positive root of
  /// k2*B^2 + k3*C*B - k1*A = 0). Used to project the DAE initial condition onto
  /// the constraint manifold (B(0)=0 is otherwise inconsistent with G=0).
  inline double ConsistentB(double k1, double k2, double k3, double a, double c)
  {
    // Positive root of k2*B^2 + k3*c*B - k1*a = 0, written in the conjugate form
    // 2*k1*a / (k3*c + sqrt(disc)) to avoid catastrophic cancellation when k3*c
    // dominates (algebraically identical to (-k3*c + sqrt(disc)) / (2*k2)).
    double disc = k3 * c * k3 * c + 4.0 * k2 * k1 * a;
    return (2.0 * k1 * a) / (k3 * c + std::sqrt(disc));
  }
}  // namespace robertson
