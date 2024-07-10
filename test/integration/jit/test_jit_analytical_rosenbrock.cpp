#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/jit/solver/jit_linear_solver.hpp>
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
  test_analytical_photolysis<builderType, stateType>(two, 1e-2);
  test_analytical_photolysis<builderType, stateType>(three, 1e-6);
  test_analytical_photolysis<builderType, stateType>(four, 1e-6);
  test_analytical_photolysis<builderType, stateType>(four_da, 1e-6);
  test_analytical_photolysis<builderType, stateType>(six_da, 1e-6);
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
  test_analytical_ternary_chemical_activation<builderType, stateType>(two, 1e-4);
  test_analytical_ternary_chemical_activation<builderType, stateType>(three, 1e-5);
  test_analytical_ternary_chemical_activation<builderType, stateType>(four, 1e-5);
  test_analytical_ternary_chemical_activation<builderType, stateType>(four_da, 1e-5);
  test_analytical_ternary_chemical_activation<builderType, stateType>(six_da, 1e-5);
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(two, 1e-2);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(three, 1e-3);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(four, 1e-3);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(four_da, 1e-3);
  test_analytical_stiff_ternary_chemical_activation<builderType, stateType>(six_da, 1e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, Tunneling)
{
  test_analytical_tunneling<builderType, stateType>(two, 1e-1);
  test_analytical_tunneling<builderType, stateType>(three, 1e-5);
  test_analytical_tunneling<builderType, stateType>(four, 1e-5);
  test_analytical_tunneling<builderType, stateType>(four_da, 1e-5);
  test_analytical_tunneling<builderType, stateType>(six_da, 1e-5);
}

TEST(AnalyticalExamplesJitRosenbrock, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling<builderType, stateType>(two, 1e-2);
  test_analytical_stiff_tunneling<builderType, stateType>(three, 1e-4);
  test_analytical_stiff_tunneling<builderType, stateType>(four, 1e-4);
  test_analytical_stiff_tunneling<builderType, stateType>(four_da, 1e-4);
  test_analytical_stiff_tunneling<builderType, stateType>(six_da, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Arrhenius)
{
  test_analytical_arrhenius<builderType, stateType>(two, 1e-4);
  test_analytical_arrhenius<builderType, stateType>(three, 1e-8);
  test_analytical_arrhenius<builderType, stateType>(four, 1e-8);
  test_analytical_arrhenius<builderType, stateType>(four_da, 1e-8);
  test_analytical_arrhenius<builderType, stateType>(six_da, 1e-8);
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
  test_analytical_branched<builderType, stateType>(two, 1e-4);
  test_analytical_branched<builderType, stateType>(three, 1e-4);
  test_analytical_branched<builderType, stateType>(four, 1e-4);
  test_analytical_branched<builderType, stateType>(four_da, 1e-4);
  test_analytical_branched<builderType, stateType>(six_da, 1e-4);
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
  auto rosenbrock_solver = [](auto params)
  {
    params.relative_tolerance_ = 1e-10;
    params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-2);
    return builderType(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_robertson<builderType, stateType>(solver, 2e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_robertson<builderType, stateType>(solver, 2e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_robertson<builderType, stateType>(solver, 2e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_robertson<builderType, stateType>(solver, 2e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_robertson<builderType, stateType>(solver, 2e-1);
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
  auto rosenbrock_solver = [](auto params)
  {
    params.relative_tolerance_ = 1e-13;
    params.absolute_tolerance_ = std::vector<double>(6, 1e-17);
    // this paper https://archimede.uniba.it/~testset/report/e5.pdf
    // says that the first variable should have a much looser tolerance than the other species
    params.absolute_tolerance_[0] = 1e-7;
    // these last two aren't actually provided values and we don't care how they behave
    params.absolute_tolerance_[4] = 1e-7;
    params.absolute_tolerance_[5] = 1e-7;
    return builderType(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_e5<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_e5<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_e5<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_e5<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_e5<builderType, stateType>(solver, 1e-3);
}

TEST(AnalyticalExamples, Oregonator)
{
  auto rosenbrock_solver = [](auto params)
  {
    // anything below 1e-6 is too strict for the Oregonator
    params.relative_tolerance_ = 1e-6;
    params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-2);
    return builderType(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_oregonator<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_oregonator<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_oregonator<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_oregonator<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_oregonator<builderType, stateType>(solver, 1e-3);
}

TEST(AnalyticalExamples, HIRES)
{
  auto rosenbrock_solver = [](auto params)
  {
    params.relative_tolerance_ = 1e-6;
    params.absolute_tolerance_ = std::vector<double>(8, params.relative_tolerance_ * 1e-2);
    return builderType(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_hires<builderType, stateType>(solver, 5e-2);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_hires<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_hires<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_hires<builderType, stateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_hires<builderType, stateType>(solver, 1e-3);
}