#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"
#include "hires.hpp"
#include "oregonator.hpp"

#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>

#include <gtest/gtest.h>

using SparseMatrixTest = micm::SparseMatrix<double>;

auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(backward_euler, 1e-4);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(backward_euler, 1e-4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(backward_euler, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(backward_euler, 1e-3);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_stiff_arrhenius(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(backward_euler, 0.05);
}

using LinearSolverTest = micm::LinearSolver<SparseMatrixTest, micm::LuDecomposition>;

template<class RatesPolicy>
using BackwardEulerTest = micm::BackwardEuler<RatesPolicy, LinearSolverTest>;

TEST(AnalyticalExamples, HIRES)
{
  using HIRESTest = HIRES<micm::Matrix<double>, SparseMatrixTest>;

  auto backward_euler_solver = HIRESTest::template CreateSolver<BackwardEulerTest<HIRESTest>, LinearSolverTest>(
      micm::BackwardEulerSolverParameters(), 1);

  test_analytical_hires(backward_euler_solver, 1e-1);
}

TEST(AnalyticalExamples, E5)
{
  test_analytical_e5(backward_euler, 1e-3);
}
