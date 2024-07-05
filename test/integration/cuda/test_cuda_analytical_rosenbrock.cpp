// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"
#include "../hires.hpp"

#include <micm/cuda/solver/cuda_solver_builder.hpp>
#include <micm/cuda/solver/cuda_solver_parameters.hpp>
#include <micm/cuda/solver/cuda_state.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <gtest/gtest.h>

constexpr std::size_t L = 1;
using builderType = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>;
using stateType = micm::CudaState<builderType::DenseMatrixPolicyType, builderType::SparseMatrixPolicyType>;

auto two = builderType(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
auto three = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four = builderType(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
auto four_da = builderType(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
auto six_da = builderType(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

auto copy_to_device = [](auto& state) -> void { state.SyncInputsToDevice(); };
auto copy_to_host = [](auto& state) -> void { state.SyncOutputsToHost(); };

// TEST(AnalyticalExamplesCudaRosenbrock, Troe)
// {
//   test_analytical_troe<builderType, stateType>(two, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_troe<builderType, stateType>(three, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_troe<builderType, stateType>(four, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_troe<builderType, stateType>(four_da, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_troe<builderType, stateType>(six_da, 1e-8, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, TroeSuperStiffButAnalytical)
// {
//   test_analytical_stiff_troe<builderType, stateType>(two, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_troe<builderType, stateType>(three, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_troe<builderType, stateType>(four, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_troe<builderType, stateType>(four_da, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_troe<builderType, stateType>(six_da, 1e-3, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, Photolysis)
// {
//   test_analytical_photolysis<builderType, stateType>(two, 1e-2, copy_to_device, copy_to_host);
//   test_analytical_photolysis<builderType, stateType>(three, 1e-6, copy_to_device, copy_to_host);
//   test_analytical_photolysis<builderType, stateType>(four, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_photolysis<builderType, stateType>(four_da, 1e-6, copy_to_device, copy_to_host);
//   test_analytical_photolysis<builderType, stateType>(six_da, 1e-8, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, PhotolysisSuperStiffButAnalytical)
// {
//   test_analytical_stiff_photolysis<builderType, stateType>(two, 1e-2, copy_to_device, copy_to_host);
//   test_analytical_stiff_photolysis<builderType, stateType>(three, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_photolysis<builderType, stateType>(four, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_stiff_photolysis<builderType, stateType>(four_da, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_stiff_photolysis<builderType, stateType>(six_da, 1e-5, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivation)
// {
//   test_analytical_ternary_chemical_activation<builderType, stateType>(two, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_ternary_chemical_activation<builderType, stateType>(three, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_ternary_chemical_activation<builderType, stateType>(four, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_ternary_chemical_activation<builderType, stateType>(four_da, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_ternary_chemical_activation<builderType, stateType>(six_da, 1e-5, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
// {
//   test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(two, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(three, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(four, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(four_da, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(six_da, 1e-3, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, Tunneling)
// {
//   test_analytical_tunneling<builderType, stateType>(two, 1e-1, copy_to_device, copy_to_host);
//   test_analytical_tunneling<builderType, stateType>(three, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_tunneling<builderType, stateType>(four, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_tunneling<builderType, stateType>(four_da, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_tunneling<builderType, stateType>(six_da, 1e-8, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, TunnelingSuperStiffButAnalytical)
// {
//   test_analytical_stiff_tunneling<builderType, stateType>(two, 1e-2, copy_to_device, copy_to_host);
//   test_analytical_stiff_tunneling<builderType, stateType>(three, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_tunneling<builderType, stateType>(four, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_tunneling<builderType, stateType>(four_da, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_tunneling<builderType, stateType>(six_da, 1e-4, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, Arrhenius)
// {
//   test_analytical_arrhenius<builderType, stateType>(two, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_arrhenius<builderType, stateType>(three, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_arrhenius<builderType, stateType>(four, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_arrhenius<builderType, stateType>(four_da, 1e-8, copy_to_device, copy_to_host);
//   test_analytical_arrhenius<builderType, stateType>(six_da, 1e-8, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, ArrheniusSuperStiffButAnalytical)
// {
//   test_analytical_stiff_arrhenius<builderType, stateType>(two, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_arrhenius<builderType, stateType>(three, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_arrhenius<builderType, stateType>(four, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_arrhenius<builderType, stateType>(four_da, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_stiff_arrhenius<builderType, stateType>(six_da, 1e-4, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, Branched)
// {
//   test_analytical_branched<builderType, stateType>(two, 1e-4, copy_to_device, copy_to_host);
//   test_analytical_branched<builderType, stateType>(three, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_branched<builderType, stateType>(four, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_branched<builderType, stateType>(four_da, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_branched<builderType, stateType>(six_da, 1e-5, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, BranchedSuperStiffButAnalytical)
// {
//   test_analytical_stiff_branched<builderType, stateType>(two, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_branched<builderType, stateType>(three, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_branched<builderType, stateType>(four, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_branched<builderType, stateType>(four_da, 1e-3, copy_to_device, copy_to_host);
//   test_analytical_stiff_branched<builderType, stateType>(six_da, 1e-3, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, Robertson)
// {
//   auto rosenbrock_solver = [](auto params) {
//     params.relative_tolerance_ = 1e-10;
//     params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-2);
//     return builderType(params);
//   };

//   auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
//   test_analytical_robertson<builderType, stateType>(solver, 2e-1, copy_to_device, copy_to_host);
  
//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
//   test_analytical_robertson<builderType, stateType>(solver, 2e-1, copy_to_device, copy_to_host);

//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
//   test_analytical_robertson<builderType, stateType>(solver, 2e-1, copy_to_device, copy_to_host);

//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
//   test_analytical_robertson<builderType, stateType>(solver, 2e-1, copy_to_device, copy_to_host);

//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
//   test_analytical_robertson<builderType, stateType>(solver, 2e-1, copy_to_device, copy_to_host);
// }

// TEST(AnalyticalExamplesCudaRosenbrock, SurfaceRxn)
// {
//   test_analytical_surface_rxn<builderType, stateType>(two, 1e-2, copy_to_device, copy_to_host);
//   test_analytical_surface_rxn<builderType, stateType>(three, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_surface_rxn<builderType, stateType>(four, 1e-6, copy_to_device, copy_to_host);
//   test_analytical_surface_rxn<builderType, stateType>(four_da, 1e-5, copy_to_device, copy_to_host);
//   test_analytical_surface_rxn<builderType, stateType>(six_da, 1e-7, copy_to_device, copy_to_host);
// }

TEST(AnalyticalExamplesCudaRosenbrock, E5)
{
  auto rosenbrock_solver = [](auto params) {
    params.relative_tolerance_ = 1e-11;
    params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-2);
    return builderType(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  // test_analytical_e5<builderType, stateType>(solver, 1e-5, copy_to_device, copy_to_host);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_e5<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);

  // solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  // test_analytical_e5<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);

  // solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  // test_analytical_e5<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);

  // solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  // test_analytical_e5<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);
}

// TEST(AnalyticalExamplesCudaRosenbrock, Oregonator)
// {
//   auto rosenbrock_solver = [](auto params) {
//     // anything below 1e-6 is too strict for the Oregonator
//     params.relative_tolerance_ = 1e-6;
//     params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-2);
//     return builderType(params);
//   };

//   auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
//   test_analytical_oregonator<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);

//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
//   test_analytical_oregonator<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);

//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
//   test_analytical_oregonator<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);

//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
//   test_analytical_oregonator<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);

//   solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
//   test_analytical_oregonator<builderType, stateType>(solver, 1e-3, copy_to_device, copy_to_host);
// }

// using LinearSolverTest = micm::CudaLinearSolver<builderType::SparseMatrixPolicyType, micm::CudaLuDecomposition>;
// template<class RatesPolicy>
// using RosenbrockTest = micm::CudaRosenbrockSolver<RatesPolicy, LinearSolverTest>;

// TEST(AnalyticalExamples, HIRES)
// {
//   using HIRESTest = HIRES<builderType::DenseMatrixPolicyType, builderType::SparseMatrixPolicyType>;

//   auto rosenbrock_solver = [](auto params) {
//     return HIRESTest::CreateSolver<RosenbrockTest<HIRESTest>, LinearSolverTest>(params, 1);
//   };

//   auto two_stage_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
//   auto three_stage_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
//   auto four_stage_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
//   auto four_stage_da_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
//   auto six_stage_da_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

//   test_analytical_hires(two_stage_solver, 1e-3);
//   test_analytical_hires(three_stage_solver, 1e-5);
//   test_analytical_hires(four_stage_solver, 1e-5);
//   test_analytical_hires(four_stage_da_solver, 1e-4);
//   test_analytical_hires(six_stage_da_solver, 1e-5);
// }