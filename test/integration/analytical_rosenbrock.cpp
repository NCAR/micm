#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"
#include "e5.hpp"
#include "hires.hpp"
#include "oregonator.hpp"

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>

#include <gtest/gtest.h>

using SparseMatrixTest = micm::SparseMatrix<double>;

auto rosenbrock_2stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
auto rosenbrock_3stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_4stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
auto rosenbrock_4stage_da = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
auto rosenbrock_6stage_da = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(rosenbrock_2stage, 1e-7);
  test_analytical_troe(rosenbrock_3stage);
  test_analytical_troe(rosenbrock_4stage);
  test_analytical_troe(rosenbrock_4stage_da);
  test_analytical_troe(rosenbrock_6stage_da);
  test_analytical_troe(backward_euler, 1e-4);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_troe(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_troe(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_troe(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_troe(rosenbrock_6stage_da, 1e-4);
  test_analytical_stiff_troe(backward_euler, 1e-4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(rosenbrock_2stage, 1e-5);
  test_analytical_photolysis(rosenbrock_3stage);
  test_analytical_photolysis(rosenbrock_4stage);
  test_analytical_photolysis(rosenbrock_4stage_da);
  test_analytical_photolysis(rosenbrock_6stage_da);
  test_analytical_photolysis(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(rosenbrock_2stage, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_photolysis(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_4stage_da, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_6stage_da, 1e-5);
  test_analytical_stiff_photolysis(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(rosenbrock_2stage);
  test_analytical_ternary_chemical_activation(rosenbrock_3stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage_da);
  test_analytical_ternary_chemical_activation(rosenbrock_6stage_da);
  test_analytical_ternary_chemical_activation(backward_euler, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_2stage, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_3stage, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage_da, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_6stage_da, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(rosenbrock_2stage, 1e-5);
  test_analytical_tunneling(rosenbrock_3stage);
  test_analytical_tunneling(rosenbrock_4stage);
  test_analytical_tunneling(rosenbrock_4stage_da);
  test_analytical_tunneling(rosenbrock_6stage_da);
  test_analytical_tunneling(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_tunneling(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_6stage_da, 1e-5);
  test_analytical_stiff_tunneling(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(rosenbrock_2stage, 1e-5);
  test_analytical_arrhenius(rosenbrock_3stage);
  test_analytical_arrhenius(rosenbrock_4stage);
  test_analytical_arrhenius(rosenbrock_4stage_da);
  test_analytical_arrhenius(rosenbrock_6stage_da);
  test_analytical_arrhenius(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-5);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-5);
  test_analytical_stiff_arrhenius(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-3);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 1e-3);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-5);
  test_analytical_stiff_arrhenius(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(rosenbrock_2stage, 1e-2);
  test_analytical_stiff_branched(rosenbrock_3stage, 1e-3);
  test_analytical_stiff_branched(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_branched(rosenbrock_4stage_da, 1e-3);
  test_analytical_stiff_branched(rosenbrock_6stage_da, 1e-5);
  test_analytical_stiff_branched(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson(rosenbrock_2stage, 1e-1);
  test_analytical_robertson(rosenbrock_3stage, 1e-1);
  test_analytical_robertson(rosenbrock_4stage, 1e-1);
  test_analytical_robertson(rosenbrock_4stage_da, 1e-1);
  test_analytical_robertson(rosenbrock_6stage_da, 1e-1);
  test_analytical_robertson(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(rosenbrock_2stage, 1e-2);
  test_analytical_surface_rxn(rosenbrock_3stage, 1e-5);
  test_analytical_surface_rxn(rosenbrock_4stage, 1e-6);
  test_analytical_surface_rxn(rosenbrock_4stage_da, 1e-5);
  test_analytical_surface_rxn(rosenbrock_6stage_da, 1e-7);
  test_analytical_surface_rxn(backward_euler, 0.05);
}

using LinearSolverTest = micm::LinearSolver<SparseMatrixTest, micm::LuDecomposition>;
template<class RatesPolicy>
using RosenbrockTest = micm::RosenbrockSolver<RatesPolicy, LinearSolverTest>;

template<class RatesPolicy>
using BackwardEulerTest = micm::BackwardEuler<RatesPolicy, LinearSolverTest>;

TEST(AnalyticalExamples, Oregonator)
{
  auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
  using OregonatorTest = Oregonator<micm::Matrix<double>, SparseMatrixTest>;

  auto rosenbrock_solver = [](auto params) {
    return OregonatorTest::template CreateSolver<RosenbrockTest<OregonatorTest>, LinearSolverTest>(params, 1);
  };

  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters()), 1e-2);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1e-2);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters()), 1e-3);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters()), 1e-2);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters()), 1e-3);
}

// TEST(AnalyticalExamples, HIRES)
// {
//   auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
//   using HIRESTest = HIRES<micm::Matrix<double>, SparseMatrixTest>;

//   auto rosenbrock_solver = HIRESTest::CreateSolver<RosenbrockTest<HIRESTest>, LinearSolverTest>(params, 1);
//   auto backward_euler_solver = HIRESTest::template CreateSolver<BackwardEulerTest<HIRESTest>, LinearSolverTest>(
//       micm::BackwardEulerSolverParameters(), 1);

//   test_analytical_hires(rosenbrock_solver, 1e-5);
//   test_analytical_hires(backward_euler_solver, 1e-1);
// }

// TEST(AnalyticalExamples, E5)
// {
//   auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
//   using E5Test = E5<micm::Matrix<double>, SparseMatrixTest>;

//   auto rosenbrock_solver = E5Test::CreateSolver<RosenbrockTest<E5Test>, LinearSolverTest>(params, 1);
//   auto backward_euler_solver =
//       E5Test::template CreateSolver<BackwardEulerTest<E5Test>, LinearSolverTest>(micm::BackwardEulerSolverParameters(), 1);

//   test_analytical_e5(rosenbrock_solver, 1e-5);
//   test_analytical_e5(backward_euler_solver, 1e-3);
// }
