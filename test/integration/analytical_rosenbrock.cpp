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

auto rosenbrock = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(rosenbrock);
  test_analytical_troe(backward_euler, 1e-4);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(rosenbrock, 1e-4);
  test_analytical_stiff_troe(backward_euler, 1e-4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(rosenbrock);
  test_analytical_photolysis(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(rosenbrock, 1e-4);
  test_analytical_stiff_photolysis(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(rosenbrock);
  test_analytical_ternary_chemical_activation(backward_euler, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(rosenbrock, 1e-4);
  test_analytical_stiff_ternary_chemical_activation(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(rosenbrock);
  test_analytical_tunneling(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(rosenbrock, 1e-4);
  test_analytical_stiff_tunneling(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(rosenbrock);
  test_analytical_arrhenius(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(rosenbrock, 1e-4);
  test_analytical_stiff_arrhenius(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched(rosenbrock, 1e-3);
  test_analytical_branched(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(rosenbrock, 1e-3);
  test_analytical_stiff_branched(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson(rosenbrock, 1e-1);
  test_analytical_robertson(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(rosenbrock, 1e-5);
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

  auto rosenbrock_solver = OregonatorTest::template CreateSolver<RosenbrockTest<OregonatorTest>, LinearSolverTest>(params, 1);
  auto backward_euler_params = micm::BackwardEulerSolverParameters();
  double abolute_tolerance = 1e-16;
  backward_euler_params.absolute_tolerance_ = { abolute_tolerance, abolute_tolerance, abolute_tolerance };
  backward_euler_params.relative_tolerance_ = 1e-10;
  auto backward_euler_solver = OregonatorTest::template CreateSolver<BackwardEulerTest<OregonatorTest>, LinearSolverTest>(backward_euler_params, 1);

  test_analytical_oregonator(rosenbrock_solver, 1e-3);
  // test_analytical_oregonator(backward_euler_solver, 1e-1);
}

TEST(AnalyticalExamples, HIRES)
{
  auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
  using HIRESTest = HIRES<micm::Matrix<double>, SparseMatrixTest>;

  auto rosenbrock_solver = HIRESTest::CreateSolver<RosenbrockTest<HIRESTest>, LinearSolverTest>(params, 1);
  auto backward_euler_solver = HIRESTest::template CreateSolver<BackwardEulerTest<HIRESTest>, LinearSolverTest>(micm::BackwardEulerSolverParameters(), 1);

  test_analytical_hires(rosenbrock_solver, 1e-5);
  test_analytical_hires(backward_euler_solver, 1e-1);
}

TEST(AnalyticalExamples, E5)
{
  auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
  using E5Test = E5<micm::Matrix<double>, SparseMatrixTest>;

  auto rosenbrock_solver = E5Test::CreateSolver<RosenbrockTest<E5Test>, LinearSolverTest>(params, 1);
  auto backward_euler_solver = E5Test::template CreateSolver<BackwardEulerTest<E5Test>, LinearSolverTest>(micm::BackwardEulerSolverParameters(), 1);

  test_analytical_e5(rosenbrock_solver, 1e-5);
  test_analytical_e5(backward_euler_solver, 1e-3);
}
