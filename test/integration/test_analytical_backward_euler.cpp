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

template<std::size_t L>
using VectorBackwardEulerDoolittle = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorStateTypeDoolittle = micm::State<
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorBackwardEulerDolittleCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorStateTypeDoolittleCSC =
    micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>, micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorBackwardEulerMozart = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateTypeMozart = micm::State<
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorBackwardEulerMozartCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>, micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateTypeMozartCSC =
    micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>, micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorBackwardEulerMozartInPlace = micm::SolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::ProcessSet,
    micm::LuDecompositionMozartInPlace,
    micm::LinearSolverInPlace<micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>, micm::LuDecompositionMozartInPlace>,
    micm::State<
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace>>;

template<std::size_t L>
using VectorStateTypeMozartInPlace = micm::State<
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace>;

auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_1 = VectorBackwardEuler<1>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_2 = VectorBackwardEuler<2>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_3 = VectorBackwardEuler<3>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());

auto backward_euler_vector_doolittle_1 = VectorBackwardEulerDoolittle<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_2 = VectorBackwardEulerDoolittle<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_3 = VectorBackwardEulerDoolittle<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_4 = VectorBackwardEulerDoolittle<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_1 = VectorBackwardEulerDolittleCSC<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_2 = VectorBackwardEulerDolittleCSC<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_3 = VectorBackwardEulerDolittleCSC<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_4 = VectorBackwardEulerDolittleCSC<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_1 = VectorBackwardEulerMozart<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_2 = VectorBackwardEulerMozart<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_3 = VectorBackwardEulerMozart<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_4 = VectorBackwardEulerMozart<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_csc_1 = VectorBackwardEulerMozartCSC<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_csc_2 = VectorBackwardEulerMozartCSC<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_csc_3 = VectorBackwardEulerMozartCSC<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_csc_4 = VectorBackwardEulerMozartCSC<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_1 = VectorBackwardEulerMozartInPlace<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_2 = VectorBackwardEulerMozartInPlace<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_3 = VectorBackwardEulerMozartInPlace<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_4 = VectorBackwardEulerMozartInPlace<4>(micm::BackwardEulerSolverParameters());

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(backward_euler, 1e-6);
  test_analytical_troe<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-6);
  test_analytical_troe<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-6);
  test_analytical_troe<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-6);
  test_analytical_troe<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-6);
  test_analytical_troe<VectorBackwardEulerDoolittle<1>, VectorStateTypeDoolittle<1>>(
      backward_euler_vector_doolittle_1, 1e-6);
  test_analytical_troe<VectorBackwardEulerDoolittle<2>, VectorStateTypeDoolittle<2>>(
      backward_euler_vector_doolittle_2, 1e-6);
  test_analytical_troe<VectorBackwardEulerDoolittle<3>, VectorStateTypeDoolittle<3>>(
      backward_euler_vector_doolittle_3, 1e-6);
  test_analytical_troe<VectorBackwardEulerDoolittle<4>, VectorStateTypeDoolittle<4>>(
      backward_euler_vector_doolittle_4, 1e-6);
  test_analytical_troe<VectorBackwardEulerDolittleCSC<1>, VectorStateTypeDoolittleCSC<1>>(backward_euler_vector_doolittle_csc_1, 1e-6);
  test_analytical_troe<VectorBackwardEulerDolittleCSC<2>, VectorStateTypeDoolittleCSC<2>>(backward_euler_vector_doolittle_csc_2, 1e-6);
  test_analytical_troe<VectorBackwardEulerDolittleCSC<3>, VectorStateTypeDoolittleCSC<3>>(backward_euler_vector_doolittle_csc_3, 1e-6);
  test_analytical_troe<VectorBackwardEulerDolittleCSC<4>, VectorStateTypeDoolittleCSC<4>>(backward_euler_vector_doolittle_csc_4, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozart<1>, VectorStateTypeMozart<1>>(backward_euler_vector_mozart_1, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozart<2>, VectorStateTypeMozart<2>>(backward_euler_vector_mozart_2, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozart<3>, VectorStateTypeMozart<3>>(backward_euler_vector_mozart_3, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozart<4>, VectorStateTypeMozart<4>>(backward_euler_vector_mozart_4, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartCSC<1>, VectorStateTypeMozartCSC<1>>(backward_euler_vector_mozart_csc_1, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartCSC<2>, VectorStateTypeMozartCSC<2>>(backward_euler_vector_mozart_csc_2, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartCSC<3>, VectorStateTypeMozartCSC<3>>(backward_euler_vector_mozart_csc_3, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartCSC<4>, VectorStateTypeMozartCSC<4>>(backward_euler_vector_mozart_csc_4, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartInPlace<1>, VectorStateTypeMozartInPlace<1>>(
      backward_euler_vector_mozart_in_place_1, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartInPlace<2>, VectorStateTypeMozartInPlace<2>>(
      backward_euler_vector_mozart_in_place_2, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartInPlace<3>, VectorStateTypeMozartInPlace<3>>(
      backward_euler_vector_mozart_in_place_3, 1e-6);
  test_analytical_troe<VectorBackwardEulerMozartInPlace<4>, VectorStateTypeMozartInPlace<4>>(
      backward_euler_vector_mozart_in_place_4, 1e-6);
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
  test_analytical_hires<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-1);
  test_analytical_hires<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-1);
  test_analytical_hires<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-1);
  test_analytical_hires<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-1);
  test_analytical_hires<VectorBackwardEulerDoolittle<1>, VectorStateTypeDoolittle<1>>(
      backward_euler_vector_doolittle_1, 1e-1);
  test_analytical_hires<VectorBackwardEulerDoolittle<2>, VectorStateTypeDoolittle<2>>(
      backward_euler_vector_doolittle_2, 1e-1);
  test_analytical_hires<VectorBackwardEulerDoolittle<3>, VectorStateTypeDoolittle<3>>(
      backward_euler_vector_doolittle_3, 1e-1);
  test_analytical_hires<VectorBackwardEulerDoolittle<4>, VectorStateTypeDoolittle<4>>(
      backward_euler_vector_doolittle_4, 1e-1);
  test_analytical_hires<VectorBackwardEulerMozart<1>, VectorStateTypeMozart<1>>(backward_euler_vector_mozart_1, 1e-1);
  test_analytical_hires<VectorBackwardEulerMozart<2>, VectorStateTypeMozart<2>>(backward_euler_vector_mozart_2, 1e-1);
  test_analytical_hires<VectorBackwardEulerMozart<3>, VectorStateTypeMozart<3>>(backward_euler_vector_mozart_3, 1e-1);
  test_analytical_hires<VectorBackwardEulerMozart<4>, VectorStateTypeMozart<4>>(backward_euler_vector_mozart_4, 1e-1);
}

TEST(AnalyticalExamples, Oregonator)
{
  test_analytical_oregonator(backward_euler, 1e-3);
  test_analytical_oregonator<VectorBackwardEuler<1>, VectorStateType<1>>(backard_euler_vector_1, 1e-3);
  test_analytical_oregonator<VectorBackwardEuler<2>, VectorStateType<2>>(backard_euler_vector_2, 1e-3);
  test_analytical_oregonator<VectorBackwardEuler<3>, VectorStateType<3>>(backard_euler_vector_3, 1e-3);
  test_analytical_oregonator<VectorBackwardEuler<4>, VectorStateType<4>>(backard_euler_vector_4, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerDoolittle<1>, VectorStateTypeDoolittle<1>>(
      backward_euler_vector_doolittle_1, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerDoolittle<2>, VectorStateTypeDoolittle<2>>(
      backward_euler_vector_doolittle_2, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerDoolittle<3>, VectorStateTypeDoolittle<3>>(
      backward_euler_vector_doolittle_3, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerDoolittle<4>, VectorStateTypeDoolittle<4>>(
      backward_euler_vector_doolittle_4, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerMozart<1>, VectorStateTypeMozart<1>>(backward_euler_vector_mozart_1, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerMozart<2>, VectorStateTypeMozart<2>>(backward_euler_vector_mozart_2, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerMozart<3>, VectorStateTypeMozart<3>>(backward_euler_vector_mozart_3, 1e-3);
  test_analytical_oregonator<VectorBackwardEulerMozart<4>, VectorStateTypeMozart<4>>(backward_euler_vector_mozart_4, 1e-3);
}
