// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/process/cuda_process.hpp>
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/process.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/branched_rate_constant.hpp>
#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/rate_constant/troe_rate_constant.hpp>
#include <micm/process/rate_constant/tunneling_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

using namespace micm;

/// @brief Test rate constant calculation on GPU with Arrhenius rate constants
///        and parameterized reactants, comparing against CPU reference
template<class DenseMatrixPolicy>
void testCudaProcessRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  ArrheniusRateConstant rc1({ .A_ = 12.2, .C_ = 300.0 });
  ArrheniusRateConstant rc2({ .A_ = 3.5e-6, .B_ = 1.2, .C_ = -200.0, .D_ = 298.0, .E_ = 1.0e-5 });
  ArrheniusRateConstant rc3({ .A_ = 1.0 });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  // CPU reference state
  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;
  }

  // Compute expected on CPU
  Process::CalculateRateConstants(processes, cpu_state);

  // Compute on GPU
  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  // Copy results back to host
  state.rate_constants_.CopyToHost();

  // Compare with relative tolerance
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with mixed rate constant types
///        (Arrhenius, Surface, UserDefined) and parameterized reactants
template<class DenseMatrixPolicy>
void testCudaProcessMixedRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  ArrheniusRateConstant rc1({ .A_ = 12.2, .C_ = 300.0 });
  SurfaceRateConstant rc2({ .label_ = "foo_surf", .phase_species_ = gas_foo });
  UserDefinedRateConstant rc3({ .label_ = "bar_user" });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;

    double radius = get_double() * 1.0e-8;
    double number_conc = get_double() * 1.0e5;
    double user_rate = get_double() * 1.0e-2;

    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.effective radius [m]"]] = radius;
    state.custom_rate_parameters_[i_cell]
                                 [state.custom_rate_parameter_map_["foo_surf.particle number concentration [# m-3]"]] =
        number_conc;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["bar_user"]] = user_rate;

    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["foo_surf.effective radius [m]"]] = radius;
    cpu_state
        .custom_rate_parameters_[i_cell]
                                [cpu_state.custom_rate_parameter_map_["foo_surf.particle number concentration [# m-3]"]] =
        number_conc;
    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["bar_user"]] = user_rate;
  }

  // Compute expected on CPU
  Process::CalculateRateConstants(processes, cpu_state);

  // Compute on GPU
  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  // Copy results back to host
  state.rate_constants_.CopyToHost();

  // Compare with relative tolerance for floating point differences
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-8)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with Troe rate constants
template<class DenseMatrixPolicy>
void testCudaProcessTroeRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  TroeRateConstant rc1({ .k0_A_ = 4.4e-32, .k0_B_ = -1.3, .kinf_A_ = 7.5e-11, .kinf_B_ = 0.2, .Fc_ = 0.6, .N_ = 1.0 });
  TroeRateConstant rc2(
      { .k0_A_ = 2.0e-30, .k0_B_ = -3.0, .k0_C_ = -394.0, .kinf_A_ = 2.5e-11, .kinf_C_ = -50.0, .Fc_ = 0.4, .N_ = 1.2 });
  TroeRateConstant rc3(TroeRateConstantParameters{});

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with Tunneling rate constants
template<class DenseMatrixPolicy>
void testCudaProcessTunnelingRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  TunnelingRateConstant rc1({ .A_ = 2.4e-12, .B_ = 167.0, .C_ = 1.65e-18 });
  TunnelingRateConstant rc2({ .A_ = 3.0e-13, .B_ = 500.0 });
  TunnelingRateConstant rc3({ .A_ = 1.0 });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with Branched rate constants (Alkoxy branch)
template<class DenseMatrixPolicy>
void testCudaProcessBranchedAlkoxyRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  BranchedRateConstant rc1(
      { .branch_ = BranchedRateConstantParameters::Branch::Alkoxy, .X_ = 2.7e-12, .Y_ = 360.0, .a0_ = 2.0e-22, .n_ = 6 });
  BranchedRateConstant rc2(
      { .branch_ = BranchedRateConstantParameters::Branch::Alkoxy, .X_ = 1.3e-11, .Y_ = 0.0, .a0_ = 1.0e-22, .n_ = 2 });
  BranchedRateConstant rc3(
      { .branch_ = BranchedRateConstantParameters::Branch::Alkoxy, .X_ = 4.0e-13, .Y_ = -200.0, .a0_ = 3.5e-22, .n_ = 8 });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with Branched rate constants (Nitrate branch)
template<class DenseMatrixPolicy>
void testCudaProcessBranchedNitrateRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  BranchedRateConstant rc1(
      { .branch_ = BranchedRateConstantParameters::Branch::Nitrate, .X_ = 2.7e-12, .Y_ = 360.0, .a0_ = 2.0e-22, .n_ = 6 });
  BranchedRateConstant rc2(
      { .branch_ = BranchedRateConstantParameters::Branch::Nitrate, .X_ = 1.3e-11, .Y_ = 0.0, .a0_ = 1.0e-22, .n_ = 2 });
  BranchedRateConstant rc3(
      { .branch_ = BranchedRateConstantParameters::Branch::Nitrate, .X_ = 4.0e-13, .Y_ = -200.0, .a0_ = 3.5e-22, .n_ = 8 });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with Ternary Chemical Activation rate constants
template<class DenseMatrixPolicy>
void testCudaProcessTernaryChemicalActivationRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  TernaryChemicalActivationRateConstant rc1(
      { .k0_A_ = 9.7e-29, .k0_B_ = -5.6, .kinf_A_ = 9.3e-12, .kinf_B_ = -1.5, .Fc_ = 0.6, .N_ = 1.0 });
  TernaryChemicalActivationRateConstant rc2(
      { .k0_A_ = 3.0e-31,
        .k0_B_ = -3.3,
        .k0_C_ = -100.0,
        .kinf_A_ = 1.5e-12,
        .kinf_C_ = -200.0,
        .Fc_ = 0.35,
        .N_ = 0.8 });
  TernaryChemicalActivationRateConstant rc3(TernaryChemicalActivationRateConstantParameters{});

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with Reversible rate constants
template<class DenseMatrixPolicy>
void testCudaProcessReversibleRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  ReversibleRateConstant rc1({ .A_ = 3.3e-10, .C_ = -1500.0, .k_r_ = 1.0e-3 });
  ReversibleRateConstant rc2({ .A_ = 1.0e-12, .k_r_ = 5.0e-6 });
  ReversibleRateConstant rc3({ .A_ = 2.5e-8, .C_ = -800.0 });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with Surface rate constants
template<class DenseMatrixPolicy>
void testCudaProcessSurfaceRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  SurfaceRateConstant rc1({ .label_ = "surf_rxn1", .phase_species_ = gas_foo, .reaction_probability_ = 0.8 });
  SurfaceRateConstant rc2({ .label_ = "surf_rxn2", .phase_species_ = gas_foo, .reaction_probability_ = 0.1 });
  SurfaceRateConstant rc3({ .label_ = "surf_rxn3", .phase_species_ = gas_foo });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;

    double radius1 = get_double() * 1.0e-8;
    double number_conc1 = get_double() * 1.0e5;
    double radius2 = get_double() * 1.0e-8;
    double number_conc2 = get_double() * 1.0e5;
    double radius3 = get_double() * 1.0e-8;
    double number_conc3 = get_double() * 1.0e5;

    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["surf_rxn1.effective radius [m]"]] = radius1;
    state.custom_rate_parameters_[i_cell]
                                 [state.custom_rate_parameter_map_["surf_rxn1.particle number concentration [# m-3]"]] =
        number_conc1;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["surf_rxn2.effective radius [m]"]] = radius2;
    state.custom_rate_parameters_[i_cell]
                                 [state.custom_rate_parameter_map_["surf_rxn2.particle number concentration [# m-3]"]] =
        number_conc2;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["surf_rxn3.effective radius [m]"]] = radius3;
    state.custom_rate_parameters_[i_cell]
                                 [state.custom_rate_parameter_map_["surf_rxn3.particle number concentration [# m-3]"]] =
        number_conc3;

    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["surf_rxn1.effective radius [m]"]] =
        radius1;
    cpu_state
        .custom_rate_parameters_[i_cell]
                                [cpu_state.custom_rate_parameter_map_["surf_rxn1.particle number concentration [# m-3]"]] =
        number_conc1;
    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["surf_rxn2.effective radius [m]"]] =
        radius2;
    cpu_state
        .custom_rate_parameters_[i_cell]
                                [cpu_state.custom_rate_parameter_map_["surf_rxn2.particle number concentration [# m-3]"]] =
        number_conc2;
    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["surf_rxn3.effective radius [m]"]] =
        radius3;
    cpu_state
        .custom_rate_parameters_[i_cell]
                                [cpu_state.custom_rate_parameter_map_["surf_rxn3.particle number concentration [# m-3]"]] =
        number_conc3;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-8)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with UserDefined rate constants
template<class DenseMatrixPolicy>
void testCudaProcessUserDefinedRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  UserDefinedRateConstant rc1({ .label_ = "user_rate1" });
  UserDefinedRateConstant rc2({ .label_ = "user_rate2", .scaling_factor_ = 2.5 });
  UserDefinedRateConstant rc3({ .label_ = "user_rate3", .scaling_factor_ = 0.1 });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;

    double user_rate1 = get_double() * 1.0e-2;
    double user_rate2 = get_double() * 1.0e-2;
    double user_rate3 = get_double() * 1.0e-2;

    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["user_rate1"]] = user_rate1;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["user_rate2"]] = user_rate2;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["user_rate3"]] = user_rate3;

    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["user_rate1"]] = user_rate1;
    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["user_rate2"]] = user_rate2;
    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["user_rate3"]] = user_rate3;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-10)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

/// @brief Test rate constant calculation on GPU with all rate constant types combined
template<class DenseMatrixPolicy>
void testCudaProcessAllRateConstants(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  ArrheniusRateConstant rc_arrhenius({ .A_ = 12.2, .C_ = 300.0 });
  TroeRateConstant rc_troe({ .k0_A_ = 4.4e-32, .k0_B_ = -1.3, .kinf_A_ = 7.5e-11, .Fc_ = 0.6 });
  TunnelingRateConstant rc_tunneling({ .A_ = 2.4e-12, .B_ = 167.0, .C_ = 1.65e-18 });
  BranchedRateConstant rc_branched_alkoxy(
      { .branch_ = BranchedRateConstantParameters::Branch::Alkoxy, .X_ = 2.7e-12, .Y_ = 360.0, .a0_ = 2.0e-22, .n_ = 6 });
  BranchedRateConstant rc_branched_nitrate(
      { .branch_ = BranchedRateConstantParameters::Branch::Nitrate, .X_ = 1.3e-11, .Y_ = 0.0, .a0_ = 1.0e-22, .n_ = 2 });
  TernaryChemicalActivationRateConstant rc_ternary(
      { .k0_A_ = 9.7e-29, .k0_B_ = -5.6, .kinf_A_ = 9.3e-12, .kinf_B_ = -1.5, .Fc_ = 0.6 });
  ReversibleRateConstant rc_reversible({ .A_ = 3.3e-10, .C_ = -1500.0, .k_r_ = 1.0e-3 });
  SurfaceRateConstant rc_surface({ .label_ = "all_surf", .phase_species_ = gas_foo, .reaction_probability_ = 0.5 });
  UserDefinedRateConstant rc_user({ .label_ = "all_user", .scaling_factor_ = 1.5 });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc_arrhenius).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc_troe).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc_tunneling).SetPhase(gas_phase).Build();
  Process r4 =
      ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc_branched_alkoxy).SetPhase(gas_phase).Build();
  Process r5 =
      ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc_branched_nitrate).SetPhase(gas_phase).Build();
  Process r6 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc_ternary).SetPhase(gas_phase).Build();
  Process r7 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc_reversible).SetPhase(gas_phase).Build();
  Process r8 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc_surface).SetPhase(gas_phase).Build();
  Process r9 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc_user).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3, r4, r5, r6, r7, r8, r9 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  State<DenseMatrixPolicy> cpu_state{ StateParameters{
                                          .number_of_rate_constants_ = processes.size(),
                                          .variable_names_ = { "foo", "bar" },
                                          .custom_rate_parameter_labels_ = param_labels,
                                      },
                                      number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    Conditions cond;
    cond.temperature_ = get_double() * 285.0;
    cond.pressure_ = get_double() * 101100.0;
    cond.air_density_ = get_double() * 10.0;
    state.conditions_[i_cell] = cond;
    cpu_state.conditions_[i_cell] = cond;

    double radius = get_double() * 1.0e-8;
    double number_conc = get_double() * 1.0e5;
    double user_rate = get_double() * 1.0e-2;

    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["all_surf.effective radius [m]"]] = radius;
    state.custom_rate_parameters_[i_cell]
                                 [state.custom_rate_parameter_map_["all_surf.particle number concentration [# m-3]"]] =
        number_conc;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["all_user"]] = user_rate;

    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["all_surf.effective radius [m]"]] =
        radius;
    cpu_state
        .custom_rate_parameters_[i_cell]
                                [cpu_state.custom_rate_parameter_map_["all_surf.particle number concentration [# m-3]"]] =
        number_conc;
    cpu_state.custom_rate_parameters_[i_cell][cpu_state.custom_rate_parameter_map_["all_user"]] = user_rate;
  }

  Process::CalculateRateConstants(processes, cpu_state);

  CudaProcess<DenseMatrixPolicy> cuda_process(processes);
  cuda_process.CalculateRateConstants(processes, state);

  state.rate_constants_.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
    {
      double expected = cpu_state.rate_constants_[i_cell][i_rxn];
      double actual = state.rate_constants_[i_cell][i_rxn];
      if (std::abs(expected) > 1.0e-30)
      {
        EXPECT_NEAR(actual / expected, 1.0, 1.0e-8)
            << "grid cell " << i_cell << "; reaction " << i_rxn << "; expected " << expected << "; actual " << actual;
      }
      else
      {
        EXPECT_NEAR(actual, expected, 1.0e-30) << "grid cell " << i_cell << "; reaction " << i_rxn;
      }
    }
  }
}

