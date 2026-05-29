// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tropospheric O3-NOx-HOx mechanism for the full-ODE vs QSSA-DAE benchmark.
// Higher-dimensional successor to the Robertson system: it carries FOUR genuine
// fast (QSSA-eligible) radicals rather than one, so the DAE removes four stiff
// modes instead of one.
//
// Timescale hierarchy (surface air, M ~ 2.5e19 cm^-3):
//   O1D  tau ~ 1 ns     (quenched by N2/O2; reacts with H2O -> 2 OH)
//   O    tau ~ 0.02 s   (O + O2 + M -> O3)
//   OH   tau ~ 1 s      (OH + CO / CH4 / NO2)
//   HO2  tau ~ 5 s      (HO2 + NO -> OH + NO2)
//   NO/NO2 ratio ~ 100 s (photostationary; NOx total is slow)
//   O3, NOx, CO, CH4, H2O2, HNO3: hours-days  (differential/slow)
//
// Photolysis rates are FIXED per run (steady-J) and supplied as custom rate
// parameters so the benchmark can sweep a solar-intensity scale factor without
// rebuilding the solver. Thermal rate constants use Arrhenius k = A*exp(C/T)
// (and (T/300)^B where B != 0), matching MICM's ArrheniusRateConstant and the
// existing Chapman regression coefficients.
//
// Reactions (rate-constant source in brackets; J* = fixed photolysis):
//   p1: NO2            -> NO + O          jNO2                       [J]
//   p2: O3             -> O1D + O2        jO1D                       [J]
//   p3: O3             -> O + O2          jO3P                       [J]
//   r1: O1D + N2       -> O + N2          2.15e-11 exp(110/T)        [JPL/Chapman]
//   r2: O1D + O2       -> O + O2          3.30e-11 exp( 55/T)        [JPL/Chapman]
//   r3: O1D + H2O      -> 2 OH            1.63e-10 exp( 60/T)        [JPL]  (HOx source)
//   r4: O + O2 + M     -> O3 + M          6.0e-34 (T/300)^-2.4 [M]   [JPL]
//   r5: O + O3         -> 2 O2            8.0e-12 exp(-2060/T)       [JPL/Chapman]
//   r6: O3 + NO        -> NO2 + O2        3.0e-12 exp(-1500/T)       [JPL]
//   r7: OH + CO (+O2)  -> HO2 + CO2       1.5e-13                    [JPL, ~P-indep here]
//   r8: OH + CH4 (+O2) -> HO2 + ...       2.45e-12 exp(-1775/T)      [JPL]
//   r9: HO2 + NO       -> OH + NO2        3.3e-12 exp( 270/T)        [JPL]  (HOx/NOx coupling)
//   r10: OH + NO2 (+M) -> HNO3            1.0e-11                    [JPL eff., termination]
//   r11: HO2 + HO2     -> H2O2 + O2       3.0e-13 exp( 460/T)        [JPL]  (HOx termination)
#pragma once

#include <micm/CPU.hpp>

#include <string>
#include <vector>

namespace tropospheric
{
  /// Default fixed (steady) photolysis rates, clear-sky surface noon (s^-1).
  /// The benchmark scales all three by a single solar factor in its sweep.
  inline constexpr double JNO2_DEFAULT = 8.0e-3;
  inline constexpr double JO1D_DEFAULT = 3.0e-5;  // O3 + hv -> O1D + O2
  inline constexpr double JO3P_DEFAULT = 4.0e-4;  // O3 + hv -> O(3P) + O2

  /// Species whose ODE row is replaced by a QSSA algebraic constraint in the
  /// DAE formulation. Ordered fast -> slow.
  inline const std::vector<std::string> QSSA_RADICALS = { "O1D", "O", "OH", "HO2" };

  /// Thermal Arrhenius parameters, the SINGLE SOURCE OF TRUTH for both the
  /// chemistry (MakeSystem) and the QSSA constraint (which recomputes the same
  /// k via micm::CalculateArrhenius at the box temperature). Naming: AR_<rxn>.
  /// k = A*exp(C/T)*(T/300)^B.
  inline constexpr micm::ArrheniusRateConstantParameters AR_O1D_N2{ .A_ = 2.15e-11, .C_ = 110 };   // O1D+N2->O+N2
  inline constexpr micm::ArrheniusRateConstantParameters AR_O1D_O2{ .A_ = 3.30e-11, .C_ = 55 };    // O1D+O2->O+O2
  inline constexpr micm::ArrheniusRateConstantParameters AR_O1D_H2O{ .A_ = 1.63e-10, .C_ = 60 };   // O1D+H2O->2OH
  inline constexpr micm::ArrheniusRateConstantParameters AR_O_O2_M{ .A_ = 6.0e-34, .B_ = -2.4 };   // O+O2+M->O3+M
  inline constexpr micm::ArrheniusRateConstantParameters AR_O_O3{ .A_ = 8.0e-12, .C_ = -2060 };    // O+O3->2O2
  inline constexpr micm::ArrheniusRateConstantParameters AR_O3_NO{ .A_ = 3.0e-12, .C_ = -1500 };   // O3+NO->NO2+O2
  inline constexpr micm::ArrheniusRateConstantParameters AR_OH_CO{ .A_ = 1.5e-13 };                // OH+CO->HO2+CO2
  inline constexpr micm::ArrheniusRateConstantParameters AR_OH_CH4{ .A_ = 2.45e-12, .C_ = -1775 }; // OH+CH4->HO2
  inline constexpr micm::ArrheniusRateConstantParameters AR_HO2_NO{ .A_ = 3.3e-12, .C_ = 270 };    // HO2+NO->OH+NO2
  inline constexpr micm::ArrheniusRateConstantParameters AR_OH_NO2{ .A_ = 1.0e-11 };               // OH+NO2->HNO3
  inline constexpr micm::ArrheniusRateConstantParameters AR_HO2_HO2{ .A_ = 3.0e-13, .C_ = 460 };   // HO2+HO2->H2O2

  struct System
  {
    micm::Phase gas_phase;
    std::vector<micm::Process> processes;
  };

  /// Build the tropospheric O3-NOx-HOx reaction system. Thermal rate constants
  /// are baked in via Arrhenius; the three photolysis rates are user-defined
  /// (labels p1/p2/p3) and set at solve time so a steady-J sweep can vary them.
  inline System MakeSystem()
  {
    // Active species.
    auto O1D = micm::Species("O1D");
    auto O = micm::Species("O");
    auto OH = micm::Species("OH");
    auto HO2 = micm::Species("HO2");
    auto O3 = micm::Species("O3");
    auto NO = micm::Species("NO");
    auto NO2 = micm::Species("NO2");
    auto CO = micm::Species("CO");
    auto CH4 = micm::Species("CH4");
    auto H2O2 = micm::Species("H2O2");
    auto HNO3 = micm::Species("HNO3");
    auto CO2 = micm::Species("CO2");
    // Large reservoirs (effectively constant over the integration; carried as
    // species so termolecular/quench reactions see them, but their fractional
    // change is negligible).
    auto O2 = micm::Species("O2");
    auto N2 = micm::Species("N2");
    auto M = micm::Species("M");
    auto H2O = micm::Species("H2O");

    micm::Phase gas_phase{ "gas",
                           std::vector<micm::PhaseSpecies>{ O1D, O, OH, HO2, O3, NO, NO2, CO, CH4, H2O2, HNO3, CO2, O2, N2,
                                                            M, H2O } };

    using micm::ArrheniusRateConstantParameters;
    using micm::ChemicalReactionBuilder;
    using micm::StoichSpecies;
    using micm::UserDefinedRateConstantParameters;

    auto arr = [&](std::vector<micm::Species> reactants,
                   std::vector<StoichSpecies> products,
                   ArrheniusRateConstantParameters p)
    {
      return ChemicalReactionBuilder()
          .SetReactants(reactants)
          .SetProducts(products)
          .SetRateConstant(p)
          .SetPhase(gas_phase)
          .Build();
    };
    auto photo = [&](std::vector<micm::Species> reactants,
                     std::vector<StoichSpecies> products,
                     const std::string& label)
    {
      return ChemicalReactionBuilder()
          .SetReactants(reactants)
          .SetProducts(products)
          .SetRateConstant(UserDefinedRateConstantParameters{ .label_ = label })
          .SetPhase(gas_phase)
          .Build();
    };

    std::vector<micm::Process> processes;

    // --- Photolysis (fixed J, set via custom rate parameters) ---
    processes.push_back(photo({ NO2 }, { { NO, 1 }, { O, 1 } }, "p1"));            // jNO2
    processes.push_back(photo({ O3 }, { { O1D, 1 }, { O2, 1 } }, "p2"));           // jO1D
    processes.push_back(photo({ O3 }, { { O, 1 }, { O2, 1 } }, "p3"));             // jO3P

    // --- O1D chemistry ---
    processes.push_back(arr({ O1D, N2 }, { { O, 1 }, { N2, 1 } }, AR_O1D_N2));
    processes.push_back(arr({ O1D, O2 }, { { O, 1 }, { O2, 1 } }, AR_O1D_O2));
    processes.push_back(arr({ O1D, H2O }, { { OH, 2 } }, AR_O1D_H2O));  // HOx source

    // --- Odd-oxygen ---
    processes.push_back(arr({ O, O2, M }, { { O3, 1 }, { M, 1 } }, AR_O_O2_M));
    processes.push_back(arr({ O, O3 }, { { O2, 2 } }, AR_O_O3));
    processes.push_back(arr({ O3, NO }, { { NO2, 1 }, { O2, 1 } }, AR_O3_NO));

    // --- HOx cycle / coupling to NOx ---
    processes.push_back(arr({ OH, CO }, { { HO2, 1 }, { CO2, 1 } }, AR_OH_CO));
    processes.push_back(arr({ OH, CH4 }, { { HO2, 1 } }, AR_OH_CH4));
    processes.push_back(arr({ HO2, NO }, { { OH, 1 }, { NO2, 1 } }, AR_HO2_NO));

    // --- Terminations ---
    processes.push_back(arr({ OH, NO2 }, { { HNO3, 1 } }, AR_OH_NO2));
    processes.push_back(arr({ HO2, HO2 }, { { H2O2, 1 }, { O2, 1 } }, AR_HO2_HO2));

    return System{ gas_phase, processes };
  }
}  // namespace tropospheric
