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
using VectorStateType = typename VectorBackwardEuler<L>::StatePolicyType;

template<std::size_t L>
using VectorBackwardEulerDoolittle = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorStateTypeDoolittle = typename VectorBackwardEulerDoolittle<L>::StatePolicyType;

template<std::size_t L>
using VectorBackwardEulerDolittleCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorStateTypeDoolittleCSC = typename VectorBackwardEulerDolittleCSC<L>::StatePolicyType;

template<std::size_t L>
using VectorBackwardEulerMozart = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateTypeMozart = typename VectorBackwardEulerMozart<L>::StatePolicyType;

template<std::size_t L>
using VectorBackwardEulerMozartCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateTypeMozartCSC = typename VectorBackwardEulerMozartCSC<L>::StatePolicyType;

template<std::size_t L>
using VectorBackwardEulerDoolittleInPlace = micm::CpuSolverBuilderInPlace<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittleInPlace>;

template<std::size_t L>
using VectorStateTypeDoolittleInPlace = typename VectorBackwardEulerDoolittleInPlace<L>::StatePolicyType;

template<std::size_t L>
using VectorBackwardEulerMozartInPlace = micm::CpuSolverBuilderInPlace<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace>;

template<std::size_t L>
using VectorStateTypeMozartInPlace = typename VectorBackwardEulerMozartInPlace<L>::StatePolicyType;

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
auto backward_euler_vector_doolittle_in_place_1 =
    VectorBackwardEulerDoolittleInPlace<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_in_place_2 =
    VectorBackwardEulerDoolittleInPlace<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_in_place_3 =
    VectorBackwardEulerDoolittleInPlace<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_in_place_4 =
    VectorBackwardEulerDoolittleInPlace<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_1 = VectorBackwardEulerMozartInPlace<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_2 = VectorBackwardEulerMozartInPlace<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_3 = VectorBackwardEulerMozartInPlace<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_4 = VectorBackwardEulerMozartInPlace<4>(micm::BackwardEulerSolverParameters());


TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(backward_euler, 1e-6);
  test_analytical_troe(backard_euler_vector_1, 1e-6);
  test_analytical_troe(backard_euler_vector_2, 1e-6);
  test_analytical_troe(backard_euler_vector_3, 1e-6);
  test_analytical_troe(backard_euler_vector_4, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_1, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_2, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_3, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_4, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_csc_1, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_csc_2, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_csc_3, 1e-6);
  test_analytical_troe(backward_euler_vector_doolittle_csc_4, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_1, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_2, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_3, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_4, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_csc_1, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_csc_2, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_csc_3, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_csc_4, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_in_place_1, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_in_place_2, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_in_place_3, 1e-6);
  test_analytical_troe(backward_euler_vector_mozart_in_place_4, 1e-6);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(backward_euler);
  test_analytical_stiff_troe(backard_euler_vector_1);
  test_analytical_stiff_troe(backard_euler_vector_2);
  test_analytical_stiff_troe(backard_euler_vector_3);
  test_analytical_stiff_troe(backard_euler_vector_4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(backward_euler, 1e-3);
  test_analytical_photolysis(backard_euler_vector_1, 1e-3);
  test_analytical_photolysis(backard_euler_vector_2, 1e-3);
  test_analytical_photolysis(backard_euler_vector_3, 1e-3);
  test_analytical_photolysis(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(backward_euler, 1e-3);
  test_analytical_stiff_photolysis(backard_euler_vector_1, 1e-3);
  test_analytical_stiff_photolysis(backard_euler_vector_2, 1e-3);
  test_analytical_stiff_photolysis(backard_euler_vector_3, 1e-3);
  test_analytical_stiff_photolysis(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(backward_euler, 1e-5);
  test_analytical_ternary_chemical_activation(backard_euler_vector_1, 1e-5);
  test_analytical_ternary_chemical_activation(backard_euler_vector_2, 1e-5);
  test_analytical_ternary_chemical_activation(backard_euler_vector_3, 1e-5);
  test_analytical_ternary_chemical_activation(backard_euler_vector_4, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(backward_euler, 1e-2);
  test_analytical_stiff_ternary_chemical_activation(backard_euler_vector_1, 1e-2);
  test_analytical_stiff_ternary_chemical_activation(backard_euler_vector_2, 1e-2);
  test_analytical_stiff_ternary_chemical_activation(backard_euler_vector_3, 1e-2);
  test_analytical_stiff_ternary_chemical_activation(backard_euler_vector_4, 1e-2);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(backward_euler, 1e-3);
  test_analytical_tunneling(backard_euler_vector_1, 1e-3);
  test_analytical_tunneling(backard_euler_vector_2, 1e-3);
  test_analytical_tunneling(backard_euler_vector_3, 1e-3);
  test_analytical_tunneling(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(backward_euler, 1e-3);
  test_analytical_stiff_tunneling(backard_euler_vector_1, 1e-3);
  test_analytical_stiff_tunneling(backard_euler_vector_2, 1e-3);
  test_analytical_stiff_tunneling(backard_euler_vector_3, 1e-3);
  test_analytical_stiff_tunneling(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(backward_euler, 1e-3);
  test_analytical_arrhenius(backard_euler_vector_1, 1e-3);
  test_analytical_arrhenius(backard_euler_vector_2, 1e-3);
  test_analytical_arrhenius(backard_euler_vector_3, 1e-3);
  test_analytical_arrhenius(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(backward_euler, 1e-3);
  test_analytical_stiff_arrhenius(backard_euler_vector_1, 1e-3);
  test_analytical_stiff_arrhenius(backard_euler_vector_2, 1e-3);
  test_analytical_stiff_arrhenius(backard_euler_vector_3, 1e-3);
  test_analytical_stiff_arrhenius(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched(backward_euler, 1e-5);
  test_analytical_branched(backard_euler_vector_1, 1e-5);
  test_analytical_branched(backard_euler_vector_2, 1e-5);
  test_analytical_branched(backard_euler_vector_3, 1e-5);
  test_analytical_branched(backard_euler_vector_4, 1e-5);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(backward_euler, 1e-2);
  test_analytical_stiff_branched(backard_euler_vector_1, 1e-2);
  test_analytical_stiff_branched(backard_euler_vector_2, 1e-2);
  test_analytical_stiff_branched(backard_euler_vector_3, 1e-2);
  test_analytical_stiff_branched(backard_euler_vector_4, 1e-2);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(backward_euler, 0.05);
}

TEST(AnalyticalExamples, HIRES)
{
  test_analytical_hires(backward_euler, 1e-1);
  test_analytical_hires(backard_euler_vector_1, 1e-1);
  test_analytical_hires(backard_euler_vector_2, 1e-1);
  test_analytical_hires(backard_euler_vector_3, 1e-1);
  test_analytical_hires(backard_euler_vector_4, 1e-1);
  test_analytical_hires(backward_euler_vector_doolittle_1, 1e-1);
  test_analytical_hires(backward_euler_vector_doolittle_2, 1e-1);
  test_analytical_hires(backward_euler_vector_doolittle_3, 1e-1);
  test_analytical_hires(backward_euler_vector_doolittle_4, 1e-1);
  test_analytical_hires(backward_euler_vector_mozart_1, 1e-1);
  test_analytical_hires(backward_euler_vector_mozart_2, 1e-1);
  test_analytical_hires(backward_euler_vector_mozart_3, 1e-1);
  test_analytical_hires(backward_euler_vector_mozart_4, 1e-1);
}

TEST(AnalyticalExamples, Oregonator)
{
  test_analytical_oregonator(backward_euler, 1e-3);
  test_analytical_oregonator(backard_euler_vector_1, 1e-3);
  test_analytical_oregonator(backard_euler_vector_2, 1e-3);
  test_analytical_oregonator(backard_euler_vector_3, 1e-3);
  test_analytical_oregonator(backard_euler_vector_4, 1e-3);
  test_analytical_oregonator(backward_euler_vector_doolittle_1, 1e-3);
  test_analytical_oregonator(backward_euler_vector_doolittle_2, 1e-3);
  test_analytical_oregonator(backward_euler_vector_doolittle_3, 1e-3);
  test_analytical_oregonator(backward_euler_vector_doolittle_4, 1e-3);
  test_analytical_oregonator(backward_euler_vector_mozart_1, 1e-3);
  test_analytical_oregonator(backward_euler_vector_mozart_2, 1e-3);
  test_analytical_oregonator(backward_euler_vector_mozart_3, 1e-3);
  test_analytical_oregonator(backward_euler_vector_mozart_4, 1e-3);
}
