#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"
#include "../oregonator.hpp"
#include "../e5.hpp"
#include "../hires.hpp"

#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/jit/solver/jit_linear_solver.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <gtest/gtest.h>

constexpr std::size_t L = 1;
using builderType = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>;
using stateType = micm::State<builderType::DenseMatrixPolicyType, builderType::SparseMatrixPolicyType>;
auto two = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto three = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four_da = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto six_da = builderType(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

TEST(AnalyticalExamplesJitRosenbrock, Troe)
{
  test_analytical_troe<builderType, stateType>(two);
  test_analytical_troe<builderType, stateType>(three);
  test_analytical_troe<builderType, stateType>(four);
  test_analytical_troe<builderType, stateType>(four_da);
  test_analytical_troe<builderType, stateType>(six_da);
}

TEST(AnalyticalExamplesJitRosenbrock, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe<builderType, stateType>(two, 1e-4);
  test_analytical_stiff_troe<builderType, stateType>(three, 1e-4);
  test_analytical_stiff_troe<builderType, stateType>(four, 1e-4);
  test_analytical_stiff_troe<builderType, stateType>(four_da, 1e-4);
  test_analytical_stiff_troe<builderType, stateType>(six_da, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Photolysis)
{
  test_analytical_photolysis<builderType, stateType>(two);
  test_analytical_photolysis<builderType, stateType>(three);
  test_analytical_photolysis<builderType, stateType>(four);
  test_analytical_photolysis<builderType, stateType>(four_da);
  test_analytical_photolysis<builderType, stateType>(six_da);
}

TEST(AnalyticalExamplesJitRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis<builderType, stateType>(two, 1e-4);
  test_analytical_stiff_photolysis<builderType, stateType>(three, 1e-4);
  test_analytical_stiff_photolysis<builderType, stateType>(four, 1e-4);
  test_analytical_stiff_photolysis<builderType, stateType>(four_da, 1e-4);
  test_analytical_stiff_photolysis<builderType, stateType>(six_da, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation<builderType, stateType>(two);
  test_analytical_ternary_chemical_activation<builderType, stateType>(three);
  test_analytical_ternary_chemical_activation<builderType, stateType>(four);
  test_analytical_ternary_chemical_activation<builderType, stateType>(four_da);
  test_analytical_ternary_chemical_activation<builderType, stateType>(six_da);
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(two, 1e-4);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(three, 1e-4);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(four, 1e-4);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(four_da, 1e-4);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(six_da, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Tunneling)
{
  test_analytical_tunneling<builderType, stateType>(two);
  test_analytical_tunneling<builderType, stateType>(three);
  test_analytical_tunneling<builderType, stateType>(four);
  test_analytical_tunneling<builderType, stateType>(four_da);
  test_analytical_tunneling<builderType, stateType>(six_da);
}

TEST(AnalyticalExamplesJitRosenbrock, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling<builderType, stateType>(two, 1e-4);
  test_analytical_stiff_tunneling<builderType, stateType>(three, 1e-4);
  test_analytical_stiff_tunneling<builderType, stateType>(four, 1e-4);
  test_analytical_stiff_tunneling<builderType, stateType>(four_da, 1e-4);
  test_analytical_stiff_tunneling<builderType, stateType>(six_da, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Arrhenius)
{
  test_analytical_arrhenius<builderType, stateType>(two);
  test_analytical_arrhenius<builderType, stateType>(three);
  test_analytical_arrhenius<builderType, stateType>(four);
  test_analytical_arrhenius<builderType, stateType>(four_da);
  test_analytical_arrhenius<builderType, stateType>(six_da);
}

TEST(AnalyticalExamplesJitRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius<builderType, stateType>(two, 1e-4);
  test_analytical_stiff_arrhenius<builderType, stateType>(three, 1e-4);
  test_analytical_stiff_arrhenius<builderType, stateType>(four, 1e-4);
  test_analytical_stiff_arrhenius<builderType, stateType>(four_da, 1e-4);
  test_analytical_stiff_arrhenius<builderType, stateType>(six_da, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Branched)
{
  test_analytical_branched<builderType, stateType>(two, 1e-3);
  test_analytical_branched<builderType, stateType>(three, 1e-3);
  test_analytical_branched<builderType, stateType>(four, 1e-3);
  test_analytical_branched<builderType, stateType>(four_da, 1e-3);
  test_analytical_branched<builderType, stateType>(six_da, 1e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched<builderType, stateType>(two, 1e-3);
  test_analytical_stiff_branched<builderType, stateType>(three, 1e-3);
  test_analytical_stiff_branched<builderType, stateType>(four, 1e-3);
  test_analytical_stiff_branched<builderType, stateType>(four_da, 1e-3);
  test_analytical_stiff_branched<builderType, stateType>(six_da, 1e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, Robertson)
{
  test_analytical_robertson<builderType, stateType>(two, 1e-1);
  test_analytical_robertson<builderType, stateType>(three, 1e-1);
  test_analytical_robertson<builderType, stateType>(four, 1e-1);
  test_analytical_robertson<builderType, stateType>(four_da, 1e-1);
  test_analytical_robertson<builderType, stateType>(six_da, 1e-1);
}

TEST(AnalyticalExamplesJitRosenbrock, SurfaceRxn)
{
  test_analytical_surface_rxn<builderType, stateType>(two, 1e-4);
  test_analytical_surface_rxn<builderType, stateType>(three, 1e-4);
  test_analytical_surface_rxn<builderType, stateType>(four, 1e-4);
  test_analytical_surface_rxn<builderType, stateType>(four_da, 1e-4);
  test_analytical_surface_rxn<builderType, stateType>(six_da, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, E5)
{
  test_analytical_e5<builderType, stateType>(two, 1e-1);
  test_analytical_e5<builderType, stateType>(three, 1e-1);
  test_analytical_e5<builderType, stateType>(four, 1e-1);
  test_analytical_e5<builderType, stateType>(four_da, 1e-1);
  test_analytical_e5<builderType, stateType>(six_da, 1e-1);
}

using LinearSolverTest = micm::JitLinearSolver<L, builderType::SparseMatrixPolicyType, micm::JitLuDecomposition<L>>;
template<class RatesPolicy>
using RosenbrockTest = micm::JitRosenbrockSolver<RatesPolicy, LinearSolverTest>;

TEST(AnalyticalExamples, Oregonator)
{
  using OregonatorTest = Oregonator<builderType::DenseMatrixPolicyType, builderType::SparseMatrixPolicyType>;

  auto rosenbrock_solver = [](auto params) {
    return OregonatorTest::template CreateSolver<RosenbrockTest<OregonatorTest>, LinearSolverTest>(params, 1);
  };

  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters()), 1e-2);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1e-2);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters()), 1e-3);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters()), 1e-2);
  test_analytical_oregonator(rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters()), 1e-3);
}

TEST(AnalyticalExamples, HIRES)
{
  using HIRESTest = HIRES<builderType::DenseMatrixPolicyType, builderType::SparseMatrixPolicyType>;

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