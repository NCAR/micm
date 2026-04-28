// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests that the CUDA rate constant kernel produces the same values as the
// CPU batch functions for each supported rate constant type.
// Each test builds a minimal solver with one reaction of the target type,
// calls UpdateStateParameters (which dispatches to GpuCalculateRateConstants),
// downloads rate_constants_ to host, and compares against the CPU calculation.

#include <micm/GPU.hpp>
#include <micm/Process.hpp>
#include <micm/Solver.hpp>
#include <micm/Util.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/util/constants.hpp>

#include <gtest/gtest.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

constexpr double TOLERANCE = 1e-10;
constexpr std::size_t L = 1;

using GpuBuilder = micm::CudaSolverBuilderInPlace<micm::CudaRosenbrockSolverParameters, L>;

const micm::Conditions kConditions{
  .temperature_ = 272.5,   // [K]
  .pressure_ = 101325.0,   // [Pa]
  .air_density_ = 2.7e19,  // [mol m-3]
};

TEST(CudaRateConstants, Arrhenius)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::ArrheniusRateConstantParameters params{ .A_ = 1.2e-10, .B_ = 0.5, .C_ = -300.0, .D_ = 273.0, .E_ = 0.0 };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  double expected = micm::CalculateArrhenius(params, kConditions.temperature_, kConditions.pressure_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, Troe)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::TroeRateConstantParameters params{ .k0_A_ = 6.0e-34,
                                           .k0_B_ = 0.0,
                                           .k0_C_ = 0.0,
                                           .kinf_A_ = 2.3e-11,
                                           .kinf_B_ = 0.0,
                                           .kinf_C_ = -200.0,
                                           .Fc_ = 0.6,
                                           .N_ = 1.0 };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  double expected = micm::CalculateTroe(params, kConditions.temperature_, kConditions.air_density_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, TernaryChemicalActivation)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::TernaryChemicalActivationRateConstantParameters params{
    .k0_A_ = 9.0e-32, .k0_B_ = -1.5, .k0_C_ = 0.0, .kinf_A_ = 3.0e-11, .kinf_B_ = 0.0, .kinf_C_ = 0.0, .Fc_ = 0.45, .N_ = 1.0
  };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  double expected = micm::CalculateTernaryChemicalActivation(params, kConditions.temperature_, kConditions.air_density_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, BranchedAlkoxy)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::BranchedRateConstantParameters params{
    .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy, .X_ = 1.2e-14, .Y_ = 0.0, .a0_ = 0.15, .n_ = 5
  };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  // Retrieve the derived params populated by BuildFrom
  const auto& store_params = solver.GetRateConstantStore().branched_[0];
  double expected = micm::CalculateBranched(store_params, kConditions.temperature_, kConditions.air_density_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, BranchedNitrate)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::BranchedRateConstantParameters params{
    .branch_ = micm::BranchedRateConstantParameters::Branch::Nitrate, .X_ = 1.2e-14, .Y_ = 0.0, .a0_ = 0.15, .n_ = 5
  };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  const auto& store_params = solver.GetRateConstantStore().branched_[0];
  double expected = micm::CalculateBranched(store_params, kConditions.temperature_, kConditions.air_density_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, Tunneling)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::TunnelingRateConstantParameters params{ .A_ = 1.0e-10, .B_ = 1200.0, .C_ = 1.7e8 };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  double expected = micm::CalculateTunneling(params, kConditions.temperature_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, TaylorSeries)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::TaylorSeriesRateConstantParameters params{};
  params.A_ = 1.5e-11;
  params.B_ = 0.0;
  params.C_ = -400.0;
  params.D_ = 300.0;
  params.E_ = 0.0;
  params.coefficients_[0] = 1.0;
  params.coefficients_[1] = 0.002;
  params.n_coefficients_ = 2;

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  double expected = micm::CalculateTaylorSeries(params, kConditions.temperature_, kConditions.pressure_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, Reversible)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::ReversibleRateConstantParameters params{ .A_ = 2.5e-12, .C_ = -180.0, .k_r_ = 3.0e-4 };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };
  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  double expected = micm::CalculateReversible(params, kConditions.temperature_);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, UserDefined)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  micm::Phase gas_phase{ "gas", { micm::PhaseSpecies(a), micm::PhaseSpecies(b) } };

  micm::UserDefinedRateConstantParameters params{ .label_ = "my_rate", .scaling_factor_ = 2.5 };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ a })
                     .SetProducts({ micm::StoichSpecies(b, 1) })
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };

  // custom_rate_parameters_ column 0 maps to label "my_rate"
  const double custom_value = 1.3e-9;
  state.custom_rate_parameters_[0][0] = custom_value;

  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  // UserDefined: rate = custom_value * scaling_factor_
  double expected = custom_value * params.scaling_factor_;
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}

TEST(CudaRateConstants, Surface)
{
  auto foo = micm::Species("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  double diffusion_coeff = 2.3e-5;
  micm::PhaseSpecies foo_gas(foo, diffusion_coeff);
  micm::Phase gas_phase{ "gas", { foo_gas } };

  micm::SurfaceRateConstantParameters params{ .label_ = "surf_rate",
                                              .phase_species_ = foo_gas,
                                              .reaction_probability_ = 0.74 };

  auto process = micm::ChemicalReactionBuilder()
                     .SetReactants({ foo })
                     .SetProducts({})
                     .SetRateConstant(params)
                     .SetPhase(gas_phase)
                     .Build();

  auto solver = GpuBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ process })
                    .Build();

  auto state = solver.GetState(L);
  state.conditions_ = { kConditions };

  // Surface reads radius at custom_param_base_index_ and num_conc at +1
  const double radius = 1.0e-7;   // [m]
  const double num_conc = 2.5e6;  // [# m-3]
  state.custom_rate_parameters_[0][0] = radius;
  state.custom_rate_parameters_[0][1] = num_conc;

  solver.UpdateStateParameters(state);
  state.rate_constants_.CopyToHost();

  // Compute expected using the same SurfaceRateConstantData the store built
  const auto& data = solver.GetRateConstantStore().surface_[0];
  double expected = micm::CalculateSurfaceOne(data, kConditions.temperature_, radius, num_conc);
  EXPECT_NEAR(state.rate_constants_[0][0], expected, TOLERANCE * std::abs(expected));
}
