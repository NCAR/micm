#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <gtest/gtest.h>

constexpr std::size_t L = 1;

TEST(AnalyticalExamplesJitRosenbrock, Troe)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_troe(
      builder, [](auto& state) {}, [](auto& state) {});
}

TEST(AnalyticalExamplesJitRosenbrock, TroeSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_troe(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Photolysis)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_photolysis(
      builder, [](auto& state) {}, [](auto& state) {});
}

TEST(AnalyticalExamplesJitRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_photolysis(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivation)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_ternary_chemical_activation(
      builder, [](auto& state) {}, [](auto& state) {});
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_ternary_chemical_activation(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Tunneling)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_tunneling(
      builder, [](auto& state) {}, [](auto& state) {});
}

TEST(AnalyticalExamplesJitRosenbrock, TunnelingSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_tunneling(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Arrhenius)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_arrhenius(
      builder, [](auto& state) {}, [](auto& state) {});
}

TEST(AnalyticalExamplesJitRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_arrhenius(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Branched)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_branched(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, BranchedSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_branched(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, Robertson)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_robertson(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-1);
}

TEST(AnalyticalExamplesJitRosenbrock, SurfaceRxn)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_surface_rxn(
      builder, [](auto& state) {}, [](auto& state) {}, 1e-4);
}
