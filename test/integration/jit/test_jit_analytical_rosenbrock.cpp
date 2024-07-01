#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
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
