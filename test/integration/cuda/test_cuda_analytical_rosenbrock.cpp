// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/GPU.hpp>

#include <gtest/gtest.h>

template<std::size_t L>
using GpuBuilder = micm::CudaSolverBuilderInPlace<micm::CudaRosenbrockSolverParameters, L>;

constexpr std::size_t L = 3;
using builderType = GpuBuilder<L>;

auto two = builderType(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
auto three = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four = builderType(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
auto four_da = builderType(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
auto six_da = builderType(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

using builderType1Cell = GpuBuilder<1>;

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
  TestAnalyticalTroe(two, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalTroe(three, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalTroe(four, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalTroe(four_da, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalTroe(six_da, 1e-10, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TroeSuperStiffButAnalytical)
{
  TestAnalyticalStiffTroe(two, 1e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffTroe(three, 1e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffTroe(four, 1e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffTroe(four_da, 1e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffTroe(six_da, 1e-5, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Photolysis)
{
  TestAnalyticalPhotolysis(two, 2e-6, copy_to_device, copy_to_host);
  TestAnalyticalPhotolysis(three, 2e-6, copy_to_device, copy_to_host);
  TestAnalyticalPhotolysis(four, 2e-8, copy_to_device, copy_to_host);
  TestAnalyticalPhotolysis(four_da, 2e-6, copy_to_device, copy_to_host);
  TestAnalyticalPhotolysis(six_da, 2e-6, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  TestAnalyticalStiffPhotolysis(two, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffPhotolysis(three, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffPhotolysis(four, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffPhotolysis(four_da, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffPhotolysis(six_da, 2e-5, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivation)
{
  TestAnalyticalTernaryChemicalActivation(two, 1e-8, copy_to_device, copy_to_host);
  TestAnalyticalTernaryChemicalActivation(three, 1e-8, copy_to_device, copy_to_host);
  TestAnalyticalTernaryChemicalActivation(four, 1e-8, copy_to_device, copy_to_host);
  TestAnalyticalTernaryChemicalActivation(four_da, 1e-8, copy_to_device, copy_to_host);
  TestAnalyticalTernaryChemicalActivation(six_da, 1e-8, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  TestAnalyticalStiffTernaryChemicalActivation(two, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffTernaryChemicalActivation(three, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffTernaryChemicalActivation(four, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffTernaryChemicalActivation(four_da, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffTernaryChemicalActivation(six_da, 2e-3, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Tunneling)
{
  TestAnalyticalTunneling(two, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalTunneling(three, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalTunneling(four, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalTunneling(four_da, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalTunneling(six_da, 1e-6, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, TunnelingSuperStiffButAnalytical)
{
  TestAnalyticalStiffTunneling(two, 1e-4, copy_to_device, copy_to_host);
  TestAnalyticalStiffTunneling(three, 1e-4, copy_to_device, copy_to_host);
  TestAnalyticalStiffTunneling(four, 1e-4, copy_to_device, copy_to_host);
  TestAnalyticalStiffTunneling(four_da, 1e-4, copy_to_device, copy_to_host);
  TestAnalyticalStiffTunneling(six_da, 1e-4, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Arrhenius)
{
  TestAnalyticalArrhenius(two, 4e-6, copy_to_device, copy_to_host);
  TestAnalyticalArrhenius(three, 1e-9, copy_to_device, copy_to_host);
  TestAnalyticalArrhenius(four, 1e-9, copy_to_device, copy_to_host);
  TestAnalyticalArrhenius(four_da, 1e-9, copy_to_device, copy_to_host);
  TestAnalyticalArrhenius(six_da, 1e-9, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  TestAnalyticalStiffArrhenius(two, 1e-4, copy_to_device, copy_to_host);
  TestAnalyticalStiffArrhenius(three, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffArrhenius(four, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffArrhenius(four_da, 2e-5, copy_to_device, copy_to_host);
  TestAnalyticalStiffArrhenius(six_da, 1e-5, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Branched)
{
  TestAnalyticalBranched(two, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalBranched(three, 1e-13, copy_to_device, copy_to_host);
  TestAnalyticalBranched(four, 1e-13, copy_to_device, copy_to_host);
  TestAnalyticalBranched(four_da, 1e-13, copy_to_device, copy_to_host);
  TestAnalyticalBranched(six_da, 1e-13, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, BranchedSuperStiffButAnalytical)
{
  TestAnalyticalStiffBranched(two, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffBranched(three, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffBranched(four, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffBranched(four_da, 2e-3, copy_to_device, copy_to_host);
  TestAnalyticalStiffBranched(six_da, 2e-3, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, SurfaceRxn)
{
  TestAnalyticalSurfaceRxn(two_1_cell, 1e-2, copy_to_device, copy_to_host);
  TestAnalyticalSurfaceRxn(three_1_cell, 1e-5, copy_to_device, copy_to_host);
  TestAnalyticalSurfaceRxn(four_1_cell, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalSurfaceRxn(four_da_1_cell, 1e-5, copy_to_device, copy_to_host);
  TestAnalyticalSurfaceRxn(six_da_1_cell, 1e-7, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Robertson)
{
  TestAnalyticalRobertson(two_1_cell, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalRobertson(three_1_cell, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalRobertson(four_1_cell, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalRobertson(four_da_1_cell, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalRobertson(six_da_1_cell, 1e-6, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, E5)
{
  TestAnalyticalE5(two_1_cell, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalE5(three_1_cell, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalE5(four_1_cell, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalE5(four_da_1_cell, 1e-10, copy_to_device, copy_to_host);
  TestAnalyticalE5(six_da_1_cell, 1e-10, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, Oregonator)
{
  TestAnalyticalOregonator(two_1_cell, 1e-2, copy_to_device, copy_to_host);
  TestAnalyticalOregonator(three_1_cell, 1e-2, copy_to_device, copy_to_host);
  TestAnalyticalOregonator(four_1_cell, 1e-2, copy_to_device, copy_to_host);
  TestAnalyticalOregonator(four_da_1_cell, 1e-2, copy_to_device, copy_to_host);
  TestAnalyticalOregonator(six_da_1_cell, 1e-2, copy_to_device, copy_to_host);
}

TEST(AnalyticalExamplesCudaRosenbrock, HIRES)
{
  TestAnalyticalHires(two_1_cell, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalHires(three_1_cell, 1e-7, copy_to_device, copy_to_host);
  TestAnalyticalHires(four_1_cell, 1e-7, copy_to_device, copy_to_host);
  TestAnalyticalHires(four_da_1_cell, 1e-6, copy_to_device, copy_to_host);
  TestAnalyticalHires(six_da_1_cell, 1e-6, copy_to_device, copy_to_host);
}
