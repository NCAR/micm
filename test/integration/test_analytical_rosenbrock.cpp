#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"
#include "hires.hpp"

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

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(rosenbrock_2stage, 1e-7);
  test_analytical_troe(rosenbrock_3stage);
  test_analytical_troe(rosenbrock_4stage);
  test_analytical_troe(rosenbrock_4stage_da);
  test_analytical_troe(rosenbrock_6stage_da);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_troe(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_troe(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_troe(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_troe(rosenbrock_6stage_da, 1e-4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(rosenbrock_2stage, 1e-5);
  test_analytical_photolysis(rosenbrock_3stage);
  test_analytical_photolysis(rosenbrock_4stage);
  test_analytical_photolysis(rosenbrock_4stage_da);
  test_analytical_photolysis(rosenbrock_6stage_da);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(rosenbrock_2stage, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_photolysis(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_4stage_da, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_6stage_da, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(rosenbrock_2stage);
  test_analytical_ternary_chemical_activation(rosenbrock_3stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage_da);
  test_analytical_ternary_chemical_activation(rosenbrock_6stage_da);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_2stage, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_3stage, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage_da, 1e-6);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_6stage_da, 1e-6);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(rosenbrock_2stage, 1e-5);
  test_analytical_tunneling(rosenbrock_3stage);
  test_analytical_tunneling(rosenbrock_4stage);
  test_analytical_tunneling(rosenbrock_4stage_da);
  test_analytical_tunneling(rosenbrock_6stage_da);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_tunneling(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_6stage_da, 1e-5);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(rosenbrock_2stage, 1e-5);
  test_analytical_arrhenius(rosenbrock_3stage);
  test_analytical_arrhenius(rosenbrock_4stage);
  test_analytical_arrhenius(rosenbrock_4stage_da);
  test_analytical_arrhenius(rosenbrock_6stage_da);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-5);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-5);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-3);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 1e-3);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-5);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(rosenbrock_2stage, 1e-2);
  test_analytical_stiff_branched(rosenbrock_3stage, 1e-3);
  test_analytical_stiff_branched(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_branched(rosenbrock_4stage_da, 1e-3);
  test_analytical_stiff_branched(rosenbrock_6stage_da, 1e-5);
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson(rosenbrock_2stage, 1e-1);
  test_analytical_robertson(rosenbrock_3stage, 1e-1);
  test_analytical_robertson(rosenbrock_4stage, 1e-1);
  test_analytical_robertson(rosenbrock_4stage_da, 1e-1);
  test_analytical_robertson(rosenbrock_6stage_da, 1e-1);
}

TEST(AnalyticalExamples, E5)
{
  auto rosenbrock_solver = [](auto params, double tolerance = 1e-8) {
    params.relative_tolerance_ = tolerance;
    params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-2);
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_e5(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_e5(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_e5(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_e5(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_e5(solver, 1e-3);
}

TEST(AnalyticalExamples, Oregonator)
{
  auto rosenbrock_solver = [](auto params) {
    // anything below 1e-6 is too strict for the Oregonator
    params.relative_tolerance_ = 1e-6;
    params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-2);
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_oregonator(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_oregonator(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_oregonator(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_oregonator(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_oregonator(solver, 1e-3);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(rosenbrock_2stage, 1e-2);
  test_analytical_surface_rxn(rosenbrock_3stage, 1e-5);
  test_analytical_surface_rxn(rosenbrock_4stage, 1e-6);
  test_analytical_surface_rxn(rosenbrock_4stage_da, 1e-5);
  test_analytical_surface_rxn(rosenbrock_6stage_da, 1e-7);
}

using LinearSolverTest = micm::LinearSolver<SparseMatrixTest, micm::LuDecomposition>;
template<class RatesPolicy>
using RosenbrockTest = micm::RosenbrockSolver<RatesPolicy, LinearSolverTest>;

TEST(AnalyticalExamples, HIRES)
{
  using HIRESTest = HIRES<micm::Matrix<double>, SparseMatrixTest>;

  auto rosenbrock_solver = [](auto params) {
    return HIRESTest::CreateSolver<RosenbrockTest<HIRESTest>, LinearSolverTest>(params, 1);
  };

  auto two_stage_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  auto three_stage_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto four_stage_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  auto four_stage_da_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  auto six_stage_da_solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

  test_analytical_hires(two_stage_solver, 1e-3);
  test_analytical_hires(three_stage_solver, 1e-5);
  test_analytical_hires(four_stage_solver, 1e-5);
  test_analytical_hires(four_stage_da_solver, 1e-4);
  test_analytical_hires(six_stage_da_solver, 1e-5);
}