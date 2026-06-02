#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

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
  TestAnalyticalTroe(backward_euler, 1e-6);
  TestAnalyticalTroe(backard_euler_vector_1, 1e-6);
  TestAnalyticalTroe(backard_euler_vector_2, 1e-6);
  TestAnalyticalTroe(backard_euler_vector_3, 1e-6);
  TestAnalyticalTroe(backard_euler_vector_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_1, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_2, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_3, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_csc_1, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_csc_2, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_csc_3, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_csc_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_1, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_2, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_3, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_csc_1, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_csc_2, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_csc_3, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_csc_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_1, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_2, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_3, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_4, 1e-6);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  TestAnalyticalStiffTroe(backward_euler);
  TestAnalyticalStiffTroe(backard_euler_vector_1);
  TestAnalyticalStiffTroe(backard_euler_vector_2);
  TestAnalyticalStiffTroe(backard_euler_vector_3);
  TestAnalyticalStiffTroe(backard_euler_vector_4);
}

TEST(AnalyticalExamples, Photolysis)
{
  TestAnalyticalPhotolysis(backward_euler, 1e-3);
  TestAnalyticalPhotolysis(backard_euler_vector_1, 1e-3);
  TestAnalyticalPhotolysis(backard_euler_vector_2, 1e-3);
  TestAnalyticalPhotolysis(backard_euler_vector_3, 1e-3);
  TestAnalyticalPhotolysis(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  TestAnalyticalStiffPhotolysis(backward_euler, 1e-3);
  TestAnalyticalStiffPhotolysis(backard_euler_vector_1, 1e-3);
  TestAnalyticalStiffPhotolysis(backard_euler_vector_2, 1e-3);
  TestAnalyticalStiffPhotolysis(backard_euler_vector_3, 1e-3);
  TestAnalyticalStiffPhotolysis(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  TestAnalyticalTernaryChemicalActivation(backward_euler, 1e-5);
  TestAnalyticalTernaryChemicalActivation(backard_euler_vector_1, 1e-5);
  TestAnalyticalTernaryChemicalActivation(backard_euler_vector_2, 1e-5);
  TestAnalyticalTernaryChemicalActivation(backard_euler_vector_3, 1e-5);
  TestAnalyticalTernaryChemicalActivation(backard_euler_vector_4, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  TestAnalyticalStiffTernaryChemicalActivation(backward_euler, 1e-2);
  TestAnalyticalStiffTernaryChemicalActivation(backard_euler_vector_1, 1e-2);
  TestAnalyticalStiffTernaryChemicalActivation(backard_euler_vector_2, 1e-2);
  TestAnalyticalStiffTernaryChemicalActivation(backard_euler_vector_3, 1e-2);
  TestAnalyticalStiffTernaryChemicalActivation(backard_euler_vector_4, 1e-2);
}

TEST(AnalyticalExamples, Tunneling)
{
  TestAnalyticalTunneling(backward_euler, 1e-3);
  TestAnalyticalTunneling(backard_euler_vector_1, 1e-3);
  TestAnalyticalTunneling(backard_euler_vector_2, 1e-3);
  TestAnalyticalTunneling(backard_euler_vector_3, 1e-3);
  TestAnalyticalTunneling(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  TestAnalyticalStiffTunneling(backward_euler, 1e-3);
  TestAnalyticalStiffTunneling(backard_euler_vector_1, 1e-3);
  TestAnalyticalStiffTunneling(backard_euler_vector_2, 1e-3);
  TestAnalyticalStiffTunneling(backard_euler_vector_3, 1e-3);
  TestAnalyticalStiffTunneling(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, Arrhenius)
{
  TestAnalyticalArrhenius(backward_euler, 1e-3);
  TestAnalyticalArrhenius(backard_euler_vector_1, 1e-3);
  TestAnalyticalArrhenius(backard_euler_vector_2, 1e-3);
  TestAnalyticalArrhenius(backard_euler_vector_3, 1e-3);
  TestAnalyticalArrhenius(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  TestAnalyticalStiffArrhenius(backward_euler, 1e-3);
  TestAnalyticalStiffArrhenius(backard_euler_vector_1, 1e-3);
  TestAnalyticalStiffArrhenius(backard_euler_vector_2, 1e-3);
  TestAnalyticalStiffArrhenius(backard_euler_vector_3, 1e-3);
  TestAnalyticalStiffArrhenius(backard_euler_vector_4, 1e-3);
}

TEST(AnalyticalExamples, Branched)
{
  TestAnalyticalBranched(backward_euler, 1e-5);
  TestAnalyticalBranched(backard_euler_vector_1, 1e-5);
  TestAnalyticalBranched(backard_euler_vector_2, 1e-5);
  TestAnalyticalBranched(backard_euler_vector_3, 1e-5);
  TestAnalyticalBranched(backard_euler_vector_4, 1e-5);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  TestAnalyticalStiffBranched(backward_euler, 1e-2);
  TestAnalyticalStiffBranched(backard_euler_vector_1, 1e-2);
  TestAnalyticalStiffBranched(backard_euler_vector_2, 1e-2);
  TestAnalyticalStiffBranched(backard_euler_vector_3, 1e-2);
  TestAnalyticalStiffBranched(backard_euler_vector_4, 1e-2);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  TestAnalyticalSurfaceRxn(backward_euler, 0.05);
}

TEST(AnalyticalExamples, HIRES)
{
  TestAnalyticalHires(backward_euler, 1e-1);
  TestAnalyticalHires(backard_euler_vector_1, 1e-1);
  TestAnalyticalHires(backard_euler_vector_2, 1e-1);
  TestAnalyticalHires(backard_euler_vector_3, 1e-1);
  TestAnalyticalHires(backard_euler_vector_4, 1e-1);
  TestAnalyticalHires(backward_euler_vector_doolittle_1, 1e-1);
  TestAnalyticalHires(backward_euler_vector_doolittle_2, 1e-1);
  TestAnalyticalHires(backward_euler_vector_doolittle_3, 1e-1);
  TestAnalyticalHires(backward_euler_vector_doolittle_4, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_1, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_2, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_3, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_4, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_1, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_2, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_3, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_4, 1e-1);
}

TEST(AnalyticalExamples, Oregonator)
{
  // The oregonator (k2=1.6e9) is an extremely stiff oscillator with a period of ~48s.
  // One large step (H=30*tau) causes Newton to converge to the wrong attractor, giving
  // ~61000x relative error. Sub-stepping at tau/1000 (~0.00016s) lets backward Euler track
  // the slow manifold with O(H) first-order accuracy (~1.2% relative error).
  TestAnalyticalOregonator(backward_euler, 0.02);
  TestAnalyticalOregonator(backard_euler_vector_1, 0.02);
  TestAnalyticalOregonator(backard_euler_vector_2, 0.02);
  TestAnalyticalOregonator(backard_euler_vector_3, 0.02);
  TestAnalyticalOregonator(backard_euler_vector_4, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_doolittle_1, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_doolittle_2, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_doolittle_3, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_doolittle_4, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_1, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_2, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_3, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_4, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_1, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_2, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_3, 0.02);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_4, 0.02);
}
