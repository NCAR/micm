#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/cuda_solver_builder.hpp>
#include <micm/solver/cuda_solver_parameters.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <gtest/gtest.h>

constexpr std::size_t L = 1;

TEST(AnalyticalExamplesCudaRosenbrock, Troe)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_troe(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, TroeSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_troe(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, Photolysis)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_photolysis(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_photolysis(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivation)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_ternary_chemical_activation(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_ternary_chemical_activation(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, Tunneling)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_tunneling(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, TunnelingSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_tunneling(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, Arrhenius)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_arrhenius(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_arrhenius(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, Branched)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_branched(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, BranchedSuperStiffButAnalytical)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_branched(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, Robertson)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_robertson(builder);
}

TEST(AnalyticalExamplesCudaRosenbrock, SurfaceRxn)
{
  auto builder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_surface_rxn(builder);
}
