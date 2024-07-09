#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>

#include <gtest/gtest.h>

using builderType = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>;
using stateType = micm::State<builderType::DenseMatrixPolicyType, builderType::SparseMatrixPolicyType>;

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
  test_analytical_troe(rosenbrock_2stage, 1e-3);
  test_analytical_troe(rosenbrock_3stage);
  test_analytical_troe(rosenbrock_4stage);
  test_analytical_troe(rosenbrock_4stage_da);
  test_analytical_troe(rosenbrock_6stage_da);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(rosenbrock_2stage, 1e-3);
  test_analytical_stiff_troe(rosenbrock_3stage, 1e-3);
  test_analytical_stiff_troe(rosenbrock_4stage, 1e-3);
  test_analytical_stiff_troe(rosenbrock_4stage_da, 1e-3);
  test_analytical_stiff_troe(rosenbrock_6stage_da, 1e-3);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(rosenbrock_2stage, 1e-2);
  test_analytical_photolysis(rosenbrock_3stage, 1e-6);
  test_analytical_photolysis(rosenbrock_4stage, 1e-8);
  test_analytical_photolysis(rosenbrock_4stage_da, 1e-6);
  test_analytical_photolysis(rosenbrock_6stage_da, 1e-8);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(rosenbrock_2stage, 1e-2);
  test_analytical_stiff_photolysis(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_photolysis(rosenbrock_4stage, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_4stage_da, 1e-5);
  test_analytical_stiff_photolysis(rosenbrock_6stage_da, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(rosenbrock_2stage, 1e-4);
  test_analytical_ternary_chemical_activation(rosenbrock_3stage, 1e-5);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage, 1e-5);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage_da, 1e-5);
  test_analytical_ternary_chemical_activation(rosenbrock_6stage_da, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_2stage, 1e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_3stage, 1e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage, 1e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage_da, 1e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_6stage_da, 1e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(rosenbrock_2stage, 1e-1);
  test_analytical_tunneling(rosenbrock_3stage, 1e-5);
  test_analytical_tunneling(rosenbrock_4stage, 1e-8);
  test_analytical_tunneling(rosenbrock_4stage_da, 1e-5);
  test_analytical_tunneling(rosenbrock_6stage_da, 1e-8);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(rosenbrock_2stage, 1e-2);
  test_analytical_stiff_tunneling(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_6stage_da, 1e-4);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(rosenbrock_2stage, 1e-4);
  test_analytical_arrhenius(rosenbrock_3stage, 1e-8);
  test_analytical_arrhenius(rosenbrock_4stage, 1e-8);
  test_analytical_arrhenius(rosenbrock_4stage_da, 1e-8);
  test_analytical_arrhenius(rosenbrock_6stage_da, 1e-8);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-4);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-4);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(rosenbrock_2stage, 1e-3);
  test_analytical_stiff_branched(rosenbrock_3stage, 1e-3);
  test_analytical_stiff_branched(rosenbrock_4stage, 1e-3);
  test_analytical_stiff_branched(rosenbrock_4stage_da, 1e-3);
  test_analytical_stiff_branched(rosenbrock_6stage_da, 1e-3);
}

TEST(AnalyticalExamples, Robertson)
{
  auto rosenbrock_solver = [](auto params) {
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

TEST(AnalyticalExamples, E5)
{
  auto rosenbrock_solver = [](auto params) {
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

TEST(AnalyticalExamples, HIRESConfig)
{
  auto rosenbrock_solver = [](auto params) {
    params.relative_tolerance_ = 1e-6;
    params.absolute_tolerance_ = std::vector<double>(8, params.relative_tolerance_ * 1e-2);
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_hires(solver, 5e-2);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_hires(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_hires(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_hires(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_hires(solver, 1e-3);
}