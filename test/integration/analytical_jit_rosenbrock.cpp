#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/jit_solver_builder.hpp>
#include <micm/solver/jit_solver_parameters.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <gtest/gtest.h>

constexpr std::size_t L = 1;

TEST(AnalyticalExamplesJitRosenbrock, Troe)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_troe(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, TroeSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_troe(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, Photolysis)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_photolysis(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_photolysis(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivation)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_ternary_chemical_activation(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_ternary_chemical_activation(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, Tunneling)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_tunneling(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, TunnelingSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_tunneling(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, Arrhenius)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_arrhenius(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_arrhenius(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, Branched)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_branched(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, BranchedSuperStiffButAnalytical)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_branched(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, Robertson)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_robertson(builder);
}

TEST(AnalyticalExamplesJitRosenbrock, SurfaceRxn)
{
  auto builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_surface_rxn(builder);
}
