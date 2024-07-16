#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

template<std::size_t L>
using VectorBackwardEuler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters, micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;
template<std::size_t L>
using VectorStateType = micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_1 = VectorBackwardEuler<1>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_2 = VectorBackwardEuler<2>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_3 = VectorBackwardEuler<3>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());


TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(backward_euler, 5e-1);
  test_analytical_troe<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_troe<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_troe<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_troe<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(backward_euler, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_stiff_troe<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(backward_euler, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_photolysis<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(backward_euler, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_stiff_photolysis<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(backward_euler, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(backward_euler, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(backward_euler, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_tunneling<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(backward_euler, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 7e-1);
  test_analytical_stiff_tunneling<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 7e-1);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(backward_euler, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_arrhenius<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(backward_euler, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched(backward_euler, 5e-1);
  test_analytical_branched<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_branched<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_branched<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_branched<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(backward_euler, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 5e-1);
  test_analytical_stiff_branched<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 5e-1);
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson(backward_euler, 1);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(backward_euler, 0.05);
}

TEST(AnalyticalExamples, HIRES)
{
  test_analytical_hires(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, E5)
{
  test_analytical_e5(backward_euler, 1);
}
