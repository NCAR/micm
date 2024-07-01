// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/cuda/solver/cuda_solver_builder.hpp>
#include <micm/cuda/solver/cuda_solver_parameters.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <gtest/gtest.h>

constexpr std::size_t L = 1;

TEST(AnalyticalExamplesCudaRosenbrock, Troe)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_troe(
      builder,
      1e-8,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, TroeSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_troe(
      builder,
      1e-4,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, Photolysis)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_photolysis(
      builder,
      1e-8,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_photolysis(
      builder,
      1e-4,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivation)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_ternary_chemical_activation(
      builder,
      1e-8,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_ternary_chemical_activation(
      builder,
      1e-4,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, Tunneling)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_tunneling(
      builder,
      1e-8
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, TunnelingSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_tunneling(
      builder,
      1e-4,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, Arrhenius)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_arrhenius(
      builder,
      1e-8,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_arrhenius(
      builder,
      1e-4,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, Branched)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_branched(
      builder,
      1e-3
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, BranchedSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_branched(
      builder,
      1e-3,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, Robertson)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_robertson(
      builder,
      1e-1,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}

TEST(AnalyticalExamplesCudaRosenbrock, SurfaceRxn)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_surface_rxn(
      builder,
      1e-5,
      [](auto& state) -> void { state.SyncInputsToDevice(); },
      [](auto& state) -> void { state.SyncOutputsToHost(); });
}
