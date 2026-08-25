#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/Kokkos.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <type_traits>

template<micm::Index L>
using VectorBackwardEuler = micm::KokkosSolverBuilder<micm::BackwardEulerSolverParameters, L>;
template<micm::Index L>
using VectorStateType = typename VectorBackwardEuler<L>::StatePolicyType;
template<micm::Index L>
using Dense = micm::KokkosDenseMatrix<micm::Real, L>;
template<micm::Index L>
using Sparse = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>;

template<micm::Index L>
using VectorBackwardEulerDoolittle = micm::
    CpuSolverBuilder<micm::BackwardEulerSolverParameters, Dense<L>, Sparse<L>, micm::LuDecompositionDoolittle<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeDoolittle = typename VectorBackwardEulerDoolittle<L>::StatePolicyType;

template<micm::Index L>
using VectorBackwardEulerDolittleCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    Dense<L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle<
        micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>>>;

template<micm::Index L>
using VectorStateTypeDoolittleCSC = typename VectorBackwardEulerDolittleCSC<L>::StatePolicyType;

template<micm::Index L>
using VectorBackwardEulerMozart =
    micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters, Dense<L>, Sparse<L>, micm::LuDecompositionMozart<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeMozart = typename VectorBackwardEulerMozart<L>::StatePolicyType;

template<micm::Index L>
using VectorBackwardEulerMozartCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    Dense<L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart<
        micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>>>;

template<micm::Index L>
using VectorStateTypeMozartCSC = typename VectorBackwardEulerMozartCSC<L>::StatePolicyType;

template<micm::Index L>
using VectorBackwardEulerDoolittleInPlace = micm::CpuSolverBuilderInPlace<
    micm::BackwardEulerSolverParameters,
    Dense<L>,
    Sparse<L>,
    micm::LuDecompositionDoolittleInPlace<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeDoolittleInPlace = typename VectorBackwardEulerDoolittleInPlace<L>::StatePolicyType;

template<micm::Index L>
using VectorBackwardEulerMozartInPlace = micm::CpuSolverBuilderInPlace<
    micm::BackwardEulerSolverParameters,
    Dense<L>,
    Sparse<L>,
    micm::LuDecompositionMozartInPlace<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeMozartInPlace = typename VectorBackwardEulerMozartInPlace<L>::StatePolicyType;

// One vector width per LU/ordering family: the other widths are covered by the Kokkos matrix
// unit tests and test_kokkos_cpu_agreement.cpp, and instantiating them all here OOMs CI.
auto backward_euler = micm::KokkosSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());

auto backward_euler_vector_doolittle_4 = VectorBackwardEulerDoolittle<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_4 = VectorBackwardEulerDolittleCSC<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_4 = VectorBackwardEulerMozart<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_csc_4 = VectorBackwardEulerMozartCSC<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_in_place_4 =
    VectorBackwardEulerDoolittleInPlace<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_4 = VectorBackwardEulerMozartInPlace<4>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_8 = VectorBackwardEulerMozartInPlace<8>(micm::BackwardEulerSolverParameters());

TEST(AnalyticalExamples, Troe)
{
  TestAnalyticalTroe(backward_euler, 1e-6);
  TestAnalyticalTroe(backard_euler_vector_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_doolittle_csc_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_csc_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_4, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_8, 1e-6);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffTroe(backward_euler);
  TestAnalyticalStiffTroe(backard_euler_vector_4);
  TestAnalyticalStiffTroe(backward_euler_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, Photolysis)
{
  TestAnalyticalPhotolysis(backward_euler, 1e-3);
  TestAnalyticalPhotolysis(backard_euler_vector_4, 1e-3);
  TestAnalyticalPhotolysis(backward_euler_vector_mozart_in_place_8, 1e-3);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffPhotolysis(backward_euler, 1e-3);
  TestAnalyticalStiffPhotolysis(backard_euler_vector_4, 1e-3);
  TestAnalyticalStiffPhotolysis(backward_euler_vector_mozart_in_place_8, 1e-3);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  TestAnalyticalTernaryChemicalActivation(backward_euler, 1e-5);
  TestAnalyticalTernaryChemicalActivation(backard_euler_vector_4, 1e-5);
  TestAnalyticalTernaryChemicalActivation(backward_euler_vector_mozart_in_place_8, 1e-5);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffTernaryChemicalActivation(backward_euler, 1e-2);
  TestAnalyticalStiffTernaryChemicalActivation(backard_euler_vector_4, 1e-2);
  TestAnalyticalStiffTernaryChemicalActivation(backward_euler_vector_mozart_in_place_8, 1e-2);
}

TEST(AnalyticalExamples, Tunneling)
{
  TestAnalyticalTunneling(backward_euler, 1e-3);
  TestAnalyticalTunneling(backard_euler_vector_4, 1e-3);
  TestAnalyticalTunneling(backward_euler_vector_mozart_in_place_8, 1e-3);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffTunneling(backward_euler, 1e-3);
  TestAnalyticalStiffTunneling(backard_euler_vector_4, 1e-3);
  TestAnalyticalStiffTunneling(backward_euler_vector_mozart_in_place_8, 1e-3);
}

TEST(AnalyticalExamples, Arrhenius)
{
  TestAnalyticalArrhenius(backward_euler, 1e-3);
  TestAnalyticalArrhenius(backard_euler_vector_4, 1e-3);
  TestAnalyticalArrhenius(backward_euler_vector_mozart_in_place_8, 1e-3);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffArrhenius(backward_euler, 1e-3);
  TestAnalyticalStiffArrhenius(backard_euler_vector_4, 1e-3);
  TestAnalyticalStiffArrhenius(backward_euler_vector_mozart_in_place_8, 1e-3);
}

TEST(AnalyticalExamples, Branched)
{
  TestAnalyticalBranched(backward_euler, 1e-5);
  TestAnalyticalBranched(backard_euler_vector_4, 1e-5);
  TestAnalyticalBranched(backward_euler_vector_mozart_in_place_8, 1e-5);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffBranched(backward_euler, 1e-2);
  TestAnalyticalStiffBranched(backard_euler_vector_4, 1e-2);
  TestAnalyticalStiffBranched(backward_euler_vector_mozart_in_place_8, 1e-2);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  TestAnalyticalSurfaceRxn(backward_euler, 0.05);
  TestAnalyticalSurfaceRxn(backward_euler_vector_mozart_in_place_8, 0.05);
}

TEST(AnalyticalExamples, HIRES)
{
  TestAnalyticalHires(backward_euler, 1e-1);
  TestAnalyticalHires(backard_euler_vector_4, 1e-1);
  TestAnalyticalHires(backward_euler_vector_doolittle_4, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_4, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_4, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_8, 1e-1);
}

TEST(AnalyticalExamples, Oregonator)
{
  // The oregonator (k2 = 1.6e9) is an extremely stiff oscillator with a period of about 48 s.
  // One large step (H = 30 * tau) makes Newton converge to the wrong attractor and gives a
  // relative error of about 61000. Sub-steps of tau/1000 let backward Euler track the slow
  // manifold with O(H) first-order accuracy. The output interval is 30 * tau, so 30000
  // sub-steps give a step size of tau/1000, and a relative error of about 3e-3.
  constexpr micm::Index kOregonatorSubsteps = 18000;
  constexpr micm::Real kOregonatorTolerance = std::is_same_v<micm::Real, double> ? 5e-3 : 1.5e-2;

  TestAnalyticalOregonator(backward_euler, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backard_euler_vector_4, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_doolittle_4, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_mozart_4, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_4, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_8, kOregonatorTolerance, kOregonatorSubsteps);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
