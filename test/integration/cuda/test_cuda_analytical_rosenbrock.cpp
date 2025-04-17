// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/GPU.hpp>

#include <gtest/gtest.h>

constexpr std::size_t L = 3;
using builderType = micm::GpuBuilder<L>;

auto two = builderType(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
auto three = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four = builderType(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
auto four_da = builderType(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
auto six_da = builderType(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

using builderType1Cell = micm::GpuBuilder<1>;

auto two_1_cell = builderType1Cell(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
auto three_1_cell = builderType1Cell(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four_1_cell = builderType1Cell(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
auto four_da_1_cell =
    builderType1Cell(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
auto six_da_1_cell = builderType1Cell(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

auto copy_to_device = [](auto& state) -> void { state.SyncInputsToDevice(); };
auto copy_to_host = [](auto& state) -> void { state.SyncOutputsToHost(); };

TEST(AnalyticalExamplesCudaRosenbrock, Troe)
{
  test_analytical_troe(two, 1e-10, copy_to_device, copy_to_host);
  test_analytical_troe(three, 1e-10, copy_to_device, copy_to_host);
  test_analytical_troe(four, 1e-10, copy_to_device, copy_to_host);
  test_analytical_troe(four_da, 1e-10, copy_to_device, copy_to_host);
  test_analytical_troe(six_da, 1e-10, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(two, 1e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_troe(three, 1e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_troe(four, 1e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_troe(four_da, 1e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_troe(six_da, 1e-5, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Photolysis)
{
  test_analytical_photolysis(two, 2e-6, copy_to_device, copy_to_host);
  test_analytical_photolysis(three, 2e-6, copy_to_device, copy_to_host);
  test_analytical_photolysis(four, 2e-8, copy_to_device, copy_to_host);
  test_analytical_photolysis(four_da, 2e-6, copy_to_device, copy_to_host);
  test_analytical_photolysis(six_da, 2e-6, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(two, 2e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_photolysis(three, 2e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_photolysis(four, 2e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_photolysis(four_da, 2e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_photolysis(six_da, 2e-5, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(two, 1e-8, copy_to_device, copy_to_host);
  test_analytical_ternary_chemical_activation(three, 1e-8, copy_to_device, copy_to_host);
  test_analytical_ternary_chemical_activation(four, 1e-8, copy_to_device, copy_to_host);
  test_analytical_ternary_chemical_activation(four_da, 1e-8, copy_to_device, copy_to_host);
  test_analytical_ternary_chemical_activation(six_da, 1e-8, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(two, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_ternary_chemical_activation(three, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_ternary_chemical_activation(four, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_ternary_chemical_activation(four_da, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_ternary_chemical_activation(six_da, 2e-3, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Tunneling)
{
  test_analytical_tunneling(two, 2e-5, copy_to_device, copy_to_host);
  test_analytical_tunneling(three, 1e-6, copy_to_device, copy_to_host);
  test_analytical_tunneling(four, 1e-6, copy_to_device, copy_to_host);
  test_analytical_tunneling(four_da, 1e-6, copy_to_device, copy_to_host);
  test_analytical_tunneling(six_da, 1e-6, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(two, 1e-4, copy_to_device, copy_to_host);
  test_analytical_stiff_tunneling(three, 1e-4, copy_to_device, copy_to_host);
  test_analytical_stiff_tunneling(four, 1e-4, copy_to_device, copy_to_host);
  test_analytical_stiff_tunneling(four_da, 1e-4, copy_to_device, copy_to_host);
  test_analytical_stiff_tunneling(six_da, 1e-4, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Arrhenius)
{
  test_analytical_arrhenius(two, 4e-6, copy_to_device, copy_to_host);
  test_analytical_arrhenius(three, 1e-9, copy_to_device, copy_to_host);
  test_analytical_arrhenius(four, 1e-9, copy_to_device, copy_to_host);
  test_analytical_arrhenius(four_da, 1e-9, copy_to_device, copy_to_host);
  test_analytical_arrhenius(six_da, 1e-9, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(two, 1e-4, copy_to_device, copy_to_host);
  test_analytical_stiff_arrhenius(three, 2e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_arrhenius(four, 2e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_arrhenius(four_da, 2e-5, copy_to_device, copy_to_host);
  test_analytical_stiff_arrhenius(six_da, 1e-5, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Branched)
{
  test_analytical_branched(two, 1e-10, copy_to_device, copy_to_host);
  test_analytical_branched(three, 1e-13, copy_to_device, copy_to_host);
  test_analytical_branched(four, 1e-13, copy_to_device, copy_to_host);
  test_analytical_branched(four_da, 1e-13, copy_to_device, copy_to_host);
  test_analytical_branched(six_da, 1e-13, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(two, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_branched(three, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_branched(four, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_branched(four_da, 2e-3, copy_to_device, copy_to_host);
  test_analytical_stiff_branched(six_da, 2e-3, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, SurfaceRxn)
{
  test_analytical_surface_rxn(two_1_cell, 1e-2, copy_to_device, copy_to_host);
  test_analytical_surface_rxn(three_1_cell, 1e-5, copy_to_device, copy_to_host);
  test_analytical_surface_rxn(four_1_cell, 1e-6, copy_to_device, copy_to_host);
  test_analytical_surface_rxn(four_da_1_cell, 1e-5, copy_to_device, copy_to_host);
  test_analytical_surface_rxn(six_da_1_cell, 1e-7, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Robertson)
{
  test_analytical_robertson(two_1_cell, 2e-1, copy_to_device, copy_to_host);
  test_analytical_robertson(three_1_cell, 2e-1, copy_to_device, copy_to_host);
  test_analytical_robertson(four_1_cell, 2e-1, copy_to_device, copy_to_host);
  test_analytical_robertson(four_da_1_cell, 2e-1, copy_to_device, copy_to_host);
  test_analytical_robertson(six_da_1_cell, 2e-1, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, E5)
{
  test_analytical_e5(two_1_cell, 1e-3, copy_to_device, copy_to_host);
  test_analytical_e5(three_1_cell, 1e-3, copy_to_device, copy_to_host);
  test_analytical_e5(four_1_cell, 1e-3, copy_to_device, copy_to_host);
  test_analytical_e5(four_da_1_cell, 1e-3, copy_to_device, copy_to_host);
  test_analytical_e5(six_da_1_cell, 1e-3, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Oregonator)
{
  test_analytical_oregonator(two_1_cell, 2e-3, copy_to_device, copy_to_host);
  test_analytical_oregonator(three_1_cell, 2e-3, copy_to_device, copy_to_host);
  test_analytical_oregonator(four_1_cell, 2e-3, copy_to_device, copy_to_host);
  test_analytical_oregonator(four_da_1_cell, 2e-3, copy_to_device, copy_to_host);
  test_analytical_oregonator(six_da_1_cell, 2e-3, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, HIRES)
{
  test_analytical_hires(two_1_cell, 1e-6, copy_to_device, copy_to_host);
  test_analytical_hires(three_1_cell, 1e-7, copy_to_device, copy_to_host);
  test_analytical_hires(four_1_cell, 1e-7, copy_to_device, copy_to_host);
  test_analytical_hires(four_da_1_cell, 1e-6, copy_to_device, copy_to_host);
  test_analytical_hires(six_da_1_cell, 1e-6, copy_to_device, copy_to_host);
}
