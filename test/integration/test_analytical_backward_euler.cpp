#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

template<std::size_t L>
using VectorBackwardEuler = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;
template<std::size_t L>
using VectorStateType =
    micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_1 = VectorBackwardEuler<1>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_2 = VectorBackwardEuler<2>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_3 = VectorBackwardEuler<3>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(backward_euler, 1e-6);
  test_analytical_troe<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-6);
  test_analytical_troe<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-6);
  test_analytical_troe<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-6);
  test_analytical_troe<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-6);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(backward_euler);
  test_analytical_stiff_troe<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1);
  test_analytical_stiff_troe<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2);
  test_analytical_stiff_troe<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3);
  test_analytical_stiff_troe<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(backward_euler, 1e-3);
  test_analytical_photolysis<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-3);
  test_analytical_photolysis<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-3);
  test_analytical_photolysis<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-3);
  test_analytical_photolysis<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(backward_euler, 1e-3);
  test_analytical_stiff_photolysis<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-3);
  test_analytical_stiff_photolysis<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-3);
  test_analytical_stiff_photolysis<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-3);
  test_analytical_stiff_photolysis<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(backward_euler, 1e-5);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-5);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-5);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-5);
  test_analytical_ternary_chemical_activation<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(backward_euler, 1e-2);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<1>, VectorStateType<1>>(
      backard_euler_vector_1, 1e-2);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<2>, VectorStateType<2>>(
      backard_euler_vector_2, 1e-2);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<3>, VectorStateType<3>>(
      backard_euler_vector_3, 1e-2);
  test_analytical_stiff_ternary_chemical_activation<VectorBackwardEuler<4>, VectorStateType<4>>(
      backard_euler_vector_4, 1e-2);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(backward_euler, 1e-3);
  test_analytical_tunneling<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-3);
  test_analytical_tunneling<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-3);
  test_analytical_tunneling<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-3);
  test_analytical_tunneling<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(backward_euler, 1e-3);
  test_analytical_stiff_tunneling<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-3);
  test_analytical_stiff_tunneling<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-3);
  test_analytical_stiff_tunneling<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-3);
  test_analytical_stiff_tunneling<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(backward_euler, 1e-3);
  test_analytical_arrhenius<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-3);
  test_analytical_arrhenius<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-3);
  test_analytical_arrhenius<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-3);
  test_analytical_arrhenius<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(backward_euler, 1e-3);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-3);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-3);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-3);
  test_analytical_stiff_arrhenius<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched(backward_euler, 1e-5);
  test_analytical_branched<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-5);
  test_analytical_branched<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-5);
  test_analytical_branched<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-5);
  test_analytical_branched<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-5);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(backward_euler, 1e-2);
  test_analytical_stiff_branched<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-2);
  test_analytical_stiff_branched<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-2);
  test_analytical_stiff_branched<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-2);
  test_analytical_stiff_branched<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-2);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(backward_euler, 0.05);
}

TEST(AnalyticalExamples, HIRES)
{
  test_analytical_hires(backward_euler, 1e-1);
}

TEST(AnalyticalExamples, Oregonator)
{
  test_analytical_oregonator(backward_euler, 1e-3);
}
