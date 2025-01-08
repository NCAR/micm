#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/jit/solver/jit_linear_solver.hpp>
#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <gtest/gtest.h>

template<std::size_t L>
using BuilderType = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>;
template<std::size_t L>
using StateType = micm::State<
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::JitLuDecomposition<L>>;

auto two = BuilderType<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto three = BuilderType<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four = BuilderType<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto four_da = BuilderType<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto six_da = BuilderType<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto param_two = micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters();
auto param_three = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
auto param_four = micm::RosenbrockSolverParameters::FourStageRosenbrockParameters();
auto param_four_da = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
auto param_six_da = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();

TEST(AnalyticalExamplesJitRosenbrock, Troe)
{
  test_analytical_troe(BuilderType<NUM_CELLS>(param_two));
  test_analytical_troe(BuilderType<NUM_CELLS>(param_three));
  test_analytical_troe(BuilderType<NUM_CELLS>(param_four));
  test_analytical_troe(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_troe(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(BuilderType<NUM_CELLS>(param_two));
  test_analytical_stiff_troe(BuilderType<NUM_CELLS>(param_three));
  test_analytical_stiff_troe(BuilderType<NUM_CELLS>(param_four));
  test_analytical_stiff_troe(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_stiff_troe(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, Photolysis)
{
  test_analytical_photolysis(BuilderType<NUM_CELLS>(param_two));
  test_analytical_photolysis(BuilderType<NUM_CELLS>(param_three));
  test_analytical_photolysis(BuilderType<NUM_CELLS>(param_four));
  test_analytical_photolysis(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_photolysis(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(BuilderType<NUM_CELLS>(param_two));
  test_analytical_stiff_photolysis(BuilderType<NUM_CELLS>(param_three));
  test_analytical_stiff_photolysis(BuilderType<NUM_CELLS>(param_four));
  test_analytical_stiff_photolysis(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_stiff_photolysis(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_two));
  test_analytical_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_three));
  test_analytical_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_four));
  test_analytical_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_two), 2e-3);
  test_analytical_stiff_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_three), 2e-3);
  test_analytical_stiff_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_four), 2e-3);
  test_analytical_stiff_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_four_da), 2e-3);
  test_analytical_stiff_ternary_chemical_activation(BuilderType<NUM_CELLS>(param_six_da), 2e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, Tunneling)
{
  test_analytical_tunneling(BuilderType<NUM_CELLS>(param_two), 2e-5);
  test_analytical_tunneling(BuilderType<NUM_CELLS>(param_three));
  test_analytical_tunneling(BuilderType<NUM_CELLS>(param_four));
  test_analytical_tunneling(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_tunneling(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(BuilderType<NUM_CELLS>(param_two), 1e-4);
  test_analytical_stiff_tunneling(BuilderType<NUM_CELLS>(param_three), 1e-4);
  test_analytical_stiff_tunneling(BuilderType<NUM_CELLS>(param_four), 1e-4);
  test_analytical_stiff_tunneling(BuilderType<NUM_CELLS>(param_four_da), 1e-4);
  test_analytical_stiff_tunneling(BuilderType<NUM_CELLS>(param_six_da), 1e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, Arrhenius)
{
  test_analytical_arrhenius(BuilderType<NUM_CELLS>(param_two), 4e-6);
  test_analytical_arrhenius(BuilderType<NUM_CELLS>(param_three));
  test_analytical_arrhenius(BuilderType<NUM_CELLS>(param_four));
  test_analytical_arrhenius(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_arrhenius(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(BuilderType<NUM_CELLS>(param_two), 1e-4);
  test_analytical_stiff_arrhenius(BuilderType<NUM_CELLS>(param_three), 2e-5);
  test_analytical_stiff_arrhenius(BuilderType<NUM_CELLS>(param_four), 2e-5);
  test_analytical_stiff_arrhenius(BuilderType<NUM_CELLS>(param_four_da), 2e-5);
  test_analytical_stiff_arrhenius(BuilderType<NUM_CELLS>(param_six_da), 1e-5);
}

TEST(AnalyticalExamplesJitRosenbrock, Branched)
{
  test_analytical_branched(BuilderType<NUM_CELLS>(param_two), 1e-10);
  test_analytical_branched(BuilderType<NUM_CELLS>(param_three));
  test_analytical_branched(BuilderType<NUM_CELLS>(param_four));
  test_analytical_branched(BuilderType<NUM_CELLS>(param_four_da));
  test_analytical_branched(BuilderType<NUM_CELLS>(param_six_da));
}

TEST(AnalyticalExamplesJitRosenbrock, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(BuilderType<NUM_CELLS>(param_two), 2e-3);
  test_analytical_stiff_branched(BuilderType<NUM_CELLS>(param_three), 2e-3);
  test_analytical_stiff_branched(BuilderType<NUM_CELLS>(param_four), 2e-3);
  test_analytical_stiff_branched(BuilderType<NUM_CELLS>(param_four_da), 2e-3);
  test_analytical_stiff_branched(BuilderType<NUM_CELLS>(param_six_da), 2e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, SurfaceRxn)
{
  test_analytical_surface_rxn(two, 1e-4);
  test_analytical_surface_rxn(three, 1e-5);
  test_analytical_surface_rxn(four, 1e-5);
  test_analytical_surface_rxn(four_da, 1e-5);
  test_analytical_surface_rxn(six_da, 1e-5);
}

TEST(AnalyticalExamplesJitRosenbrock, Robertson)
{
  test_analytical_robertson(BuilderType<1>(param_two), 1e-1);
  test_analytical_robertson(BuilderType<1>(param_three), 1e-1);
  test_analytical_robertson(BuilderType<1>(param_four), 1e-1);
  test_analytical_robertson(BuilderType<1>(param_four_da), 1e-1);
  test_analytical_robertson(BuilderType<1>(param_six_da), 1e-1);
}

TEST(AnalyticalExamplesJitRosenbrock, E5)
{
  test_analytical_e5(BuilderType<1>(param_two), 1e-3);
  test_analytical_e5(BuilderType<1>(param_three), 1e-3);
  test_analytical_e5(BuilderType<1>(param_four), 1e-3);
  test_analytical_e5(BuilderType<1>(param_four_da), 1e-3);
  test_analytical_e5(BuilderType<1>(param_six_da), 1e-3);
}

TEST(AnalyticalExamplesJitRosenbrock, Oregonator)
{
  test_analytical_oregonator(BuilderType<1>(param_two), 4e-4);
  test_analytical_oregonator(BuilderType<1>(param_three), 4e-4);
  test_analytical_oregonator(BuilderType<1>(param_four), 4e-4);
  test_analytical_oregonator(BuilderType<1>(param_four_da), 4e-4);
  test_analytical_oregonator(BuilderType<1>(param_six_da), 4e-4);
}

TEST(AnalyticalExamplesJitRosenbrock, HIRES)
{
  test_analytical_hires(BuilderType<1>(param_two), 1e-6);
  test_analytical_hires(BuilderType<1>(param_three), 1e-7);
  test_analytical_hires(BuilderType<1>(param_four), 1e-7);
  test_analytical_hires(BuilderType<1>(param_four_da), 1e-6);
  test_analytical_hires(BuilderType<1>(param_six_da), 1e-6);
}
