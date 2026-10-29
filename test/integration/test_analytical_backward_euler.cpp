#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <type_traits>

template<micm::Index L>
using VectorBackwardEuler = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;
template<micm::Index L>
using VectorStateType = typename VectorBackwardEuler<L>::StatePolicyType;

template<micm::Index L>
using VectorBackwardEulerMozartInPlace = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace<micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>>;

template<micm::Index L>
using VectorStateTypeMozartInPlace = typename VectorBackwardEulerMozartInPlace<L>::StatePolicyType;

auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_1 = VectorBackwardEuler<1>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_2 = VectorBackwardEuler<2>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_3 = VectorBackwardEuler<3>(micm::BackwardEulerSolverParameters());
auto backard_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());

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
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_1, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_2, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_3, 1e-6);
  TestAnalyticalTroe(backward_euler_vector_mozart_in_place_4, 1e-6);
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
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

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
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

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
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

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
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

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
  // The fast equilibrium in these systems (k ~ 4e10 against h = 1) makes the 2x2 block of the
  // backward-Euler matrix exactly singular in single precision: (1+k3)(1+k4) and k3*k4 round to
  // the same float, so the pivot is 0 and the solve returns inf.
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

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
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_1, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_2, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_3, 1e-1);
  TestAnalyticalHires(backward_euler_vector_mozart_in_place_4, 1e-1);
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
  TestAnalyticalOregonator(backard_euler_vector_1, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backard_euler_vector_2, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backard_euler_vector_3, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backard_euler_vector_4, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_1, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_2, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_3, kOregonatorTolerance, kOregonatorSubsteps);
  TestAnalyticalOregonator(backward_euler_vector_mozart_in_place_4, kOregonatorTolerance, kOregonatorSubsteps);
}
