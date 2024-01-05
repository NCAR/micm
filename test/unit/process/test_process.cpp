// Copyright (C) 2023-2024 National Center for Atmospheric Research,
// SPDX-License-Identifier: Apache-2.0
#include <gtest/gtest.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>

template<template<class> class MatrixPolicy>
void testProcessUpdateState(const std::size_t number_of_grid_cells)
{
  micm::Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  micm::Species bar("bar");
  bar.parameterize_ = [](const micm::Conditions& c) { return c.air_density_ * 0.82; };

  micm::ArrheniusRateConstant rc1({ .A_ = 12.2, .C_ = 300.0 });
  micm::SurfaceRateConstant rc2({ .label_ = "foo_surf", .species_ = foo });
  micm::UserDefinedRateConstant rc3({ .label_ = "bar_user" });

  micm::Process r1 = micm::Process::create().rate_constant(rc1).reactants({ foo, bar });
  micm::Process r2 = micm::Process::create().rate_constant(rc2);
  micm::Process r3 = micm::Process::create().rate_constant(rc3);
  std::vector<micm::Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
    for (auto& label : process.rate_constant_->CustomParameters())
      param_labels.push_back(label);
  micm::State<MatrixPolicy> state{ micm::StateParameters{
      .number_of_grid_cells_ = number_of_grid_cells,
      .number_of_rate_constants_ = processes.size(),
      .variable_names_ = { "foo", "bar'" },
      .custom_rate_parameter_labels_ = param_labels,
  } };

  MatrixPolicy<double> expected_rate_constants(number_of_grid_cells, 3, 0.0);
  std::vector<double> params = { 0.0, 0.0, 0.0 };
  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    state.conditions_[i_cell].temperature_ = get_double() * 285.0;
    state.conditions_[i_cell].pressure_ = get_double() * 101100.0;
    state.conditions_[i_cell].air_density_ = get_double() * 10.0;
    params[0] = get_double() * 1.0e-8;
    params[1] = get_double() * 1.0e5;
    params[2] = get_double() * 1.0e-2;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.effective radius [m]"]] = params[0];
    state.custom_rate_parameters_[i_cell]
                                 [state.custom_rate_parameter_map_["foo_surf.particle number concentration [# m-3]"]] =
        params[1];
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["bar_user"]] = params[2];
    std::vector<double>::const_iterator param_iter = params.begin();
    expected_rate_constants[i_cell][0] =
        rc1.calculate(state.conditions_[i_cell], param_iter) * (state.conditions_[i_cell].air_density_ * 0.82);
    param_iter += rc1.SizeCustomParameters();
    expected_rate_constants[i_cell][1] = rc2.calculate(state.conditions_[i_cell], param_iter);
    param_iter += rc2.SizeCustomParameters();
    expected_rate_constants[i_cell][2] = rc3.calculate(state.conditions_[i_cell], param_iter);
    param_iter += rc3.SizeCustomParameters();
  }

  micm::Process::UpdateState(processes, state);

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
      EXPECT_EQ(state.rate_constants_[i_cell][i_rxn], expected_rate_constants[i_cell][i_rxn])
          << "grid cell " << i_cell << "; reaction " << i_rxn;
}

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

TEST(Process, Matrix)
{
  testProcessUpdateState<micm::Matrix>(5);
}

TEST(Process, VectorMatrix)
{
  testProcessUpdateState<Group1VectorMatrix>(5);
  testProcessUpdateState<Group2VectorMatrix>(5);
  testProcessUpdateState<Group3VectorMatrix>(5);
  testProcessUpdateState<Group4VectorMatrix>(5);
}

TEST(Process, SurfaceRateConstantOnlyHasOneReactant)
{
  micm::Species c("c", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  micm::Species e("e");

  micm::Phase gas_phase({ c, e });
  EXPECT_ANY_THROW(micm::Process r = micm::Process::create()
                                         .reactants({ c, c })
                                         .products({ yields(e, 1) })
                                         .rate_constant(micm::SurfaceRateConstant(
                                             { .label_ = "c", .species_ = c, .reaction_probability_ = 0.90 }))
                                         .phase(gas_phase););
}
