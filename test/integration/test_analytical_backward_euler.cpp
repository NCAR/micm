#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>


using BuilderType = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>;
using StateType = micm::State<micm::BackwardEulerTemporaryVariables<BuilderType::DenseMatrixPolicyType>, BuilderType::DenseMatrixPolicyType, BuilderType::SparseMatrixPolicyType>;

template<std::size_t L>
using VectorBackwardEuler = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;
template<std::size_t L>
using VectorStateType =
    micm::State<micm::BackwardEulerTemporaryVariables<micm::VectorMatrix<double, L>>, micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

auto backward_euler = BuilderType(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_1 = VectorBackwardEuler<1>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_2 = VectorBackwardEuler<2>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_3 = VectorBackwardEuler<3>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_troe<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_troe<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_troe<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_troe<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis<BuilderType, StateType>(backward_euler, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis<BuilderType, StateType>(backward_euler, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<1>, VectorStateType<1>>(
      backard_euler_vector_1, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<2>, VectorStateType<2>>(
      backard_euler_vector_2, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<3>, VectorStateType<3>>(
      backard_euler_vector_3, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<4>, VectorStateType<4>>(
      backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling<BuilderType, StateType>(backward_euler, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling<BuilderType, StateType>(backward_euler, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_branched<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_branched<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_branched<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_branched<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched<BuilderType, StateType>(backward_euler, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson<BuilderType, StateType>(backward_euler, 1);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn<BuilderType, StateType>(backward_euler, 0.05);
}

TEST(AnalyticalExamples, HIRES)
{
  test_analytical_hires<BuilderType, StateType>(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, E5)
{
  test_analytical_e5<BuilderType, StateType>(backward_euler, 1);
}
