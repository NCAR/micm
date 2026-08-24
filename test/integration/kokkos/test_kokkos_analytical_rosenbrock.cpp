#include "../analytical_policy.hpp"
#include "../analytical_surface_rxn_policy.hpp"

#include <micm/Kokkos.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <type_traits>

using BuilderType = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>;
using StateType = micm::State<BuilderType::DenseMatrixPolicyType, BuilderType::SparseMatrixPolicyType>;
template<micm::Index L>
using Dense = micm::KokkosDenseMatrix<micm::Real, L>;
template<micm::Index L>
using Sparse = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>;

template<micm::Index L>
using VectorRosenbrock = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, Dense<L>, Sparse<L>>;

template<micm::Index L>
using VectorRosenbrockDoolittle =
    micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, Dense<L>, Sparse<L>, micm::LuDecompositionDoolittle<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeDoolittle = micm::State<Dense<L>, Sparse<L>, micm::LuDecompositionDoolittle<Sparse<L>>>;

template<micm::Index L>
using VectorRosenbrockMozart =
    micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, Dense<L>, Sparse<L>, micm::LuDecompositionMozart<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeMozart = micm::State<Dense<L>, Sparse<L>, micm::LuDecompositionMozart<Sparse<L>>>;

template<micm::Index L>
using VectorStateType = micm::State<Dense<L>, Sparse<L>>;

template<micm::Index L>
using VectorRosenbrockDolittleCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    Dense<L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle<
        micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>>>;

template<micm::Index L>
using VectorStateTypeDoolittleCSC = typename VectorRosenbrockDolittleCSC<L>::StatePolicyType;

template<micm::Index L>
using VectorRosenbrockMozartCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    Dense<L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart<
        micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>>>;

template<micm::Index L>
using VectorStateTypeMozartCSC = typename VectorRosenbrockMozartCSC<L>::StatePolicyType;

template<micm::Index L>
using VectorRosenbrockDoolittleInPlace = micm::CpuSolverBuilderInPlace<
    micm::RosenbrockSolverParameters,
    Dense<L>,
    Sparse<L>,
    micm::LuDecompositionDoolittleInPlace<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeDoolittleInPlace = typename VectorRosenbrockDoolittleInPlace<L>::StatePolicyType;

template<micm::Index L>
using VectorRosenbrockMozartInPlace = micm::CpuSolverBuilderInPlace<
    micm::RosenbrockSolverParameters,
    Dense<L>,
    Sparse<L>,
    micm::LuDecompositionMozartInPlace<Sparse<L>>>;

template<micm::Index L>
using VectorStateTypeMozartInPlace = typename VectorRosenbrockMozartInPlace<L>::StatePolicyType;

// One vector width per LU/ordering family: the other widths are covered by the Kokkos matrix
// unit tests and test_kokkos_cpu_agreement.cpp, and instantiating them all here OOMs CI.
auto rosenbrock_2stage = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
auto rosenbrock_3stage = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_4stage = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
auto rosenbrock_4stage_da = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
auto rosenbrock_6stage_da = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

auto rosenbrock_vector_4 = VectorRosenbrock<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_doolittle_4 =
    VectorRosenbrockDoolittle<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_4 =
    VectorRosenbrockMozart<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_doolittle_csc_4 =
    VectorRosenbrockDolittleCSC<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_mozart_csc_4 =
    VectorRosenbrockMozartCSC<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_doolittle_in_place_4 =
    VectorRosenbrockDoolittleInPlace<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_mozart_in_place_4 =
    VectorRosenbrockMozartInPlace<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_in_place_8 =
    VectorRosenbrockMozartInPlace<8>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

TEST(AnalyticalExamples, Troe)
{
  TestAnalyticalTroe(rosenbrock_2stage);
  TestAnalyticalTroe(rosenbrock_3stage);
  TestAnalyticalTroe(rosenbrock_4stage);
  TestAnalyticalTroe(rosenbrock_4stage_da);
  TestAnalyticalTroe(rosenbrock_6stage_da);
  TestAnalyticalTroe(rosenbrock_vector_4);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_4);
  TestAnalyticalTroe(rosenbrock_vector_mozart_4);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalTroe(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalTroe(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalTroe(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffTroe(rosenbrock_2stage);
  TestAnalyticalStiffTroe(rosenbrock_3stage);
  TestAnalyticalStiffTroe(rosenbrock_4stage);
  TestAnalyticalStiffTroe(rosenbrock_4stage_da);
  TestAnalyticalStiffTroe(rosenbrock_6stage_da);
  TestAnalyticalStiffTroe(rosenbrock_vector_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, Photolysis)
{
  TestAnalyticalPhotolysis(rosenbrock_2stage);
  TestAnalyticalPhotolysis(rosenbrock_3stage);
  TestAnalyticalPhotolysis(rosenbrock_4stage);
  TestAnalyticalPhotolysis(rosenbrock_4stage_da);
  TestAnalyticalPhotolysis(rosenbrock_6stage_da);
  TestAnalyticalPhotolysis(rosenbrock_vector_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffPhotolysis(rosenbrock_2stage);
  TestAnalyticalStiffPhotolysis(rosenbrock_3stage);
  TestAnalyticalStiffPhotolysis(rosenbrock_4stage);
  TestAnalyticalStiffPhotolysis(rosenbrock_4stage_da);
  TestAnalyticalStiffPhotolysis(rosenbrock_6stage_da);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  TestAnalyticalTernaryChemicalActivation(rosenbrock_2stage);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_3stage);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_4stage);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_4stage_da);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_6stage_da);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_2stage, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_3stage, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_4stage, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_4stage_da, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_6stage_da, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_csc_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_8, 2e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  TestAnalyticalTunneling(rosenbrock_2stage, 2e-5);
  TestAnalyticalTunneling(rosenbrock_3stage);
  TestAnalyticalTunneling(rosenbrock_4stage);
  TestAnalyticalTunneling(rosenbrock_4stage_da);
  TestAnalyticalTunneling(rosenbrock_6stage_da);
  TestAnalyticalTunneling(rosenbrock_vector_4);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffTunneling(rosenbrock_2stage, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_3stage, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_4stage, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_4stage_da, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_6stage_da, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_csc_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_csc_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_in_place_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_in_place_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_in_place_8, 1e-4);
}

TEST(AnalyticalExamples, Arrhenius)
{
  TestAnalyticalArrhenius(rosenbrock_2stage, 4e-6);
  TestAnalyticalArrhenius(rosenbrock_3stage);
  TestAnalyticalArrhenius(rosenbrock_4stage);
  TestAnalyticalArrhenius(rosenbrock_4stage_da);
  TestAnalyticalArrhenius(rosenbrock_6stage_da);
  TestAnalyticalArrhenius(rosenbrock_vector_4);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffArrhenius(rosenbrock_2stage, 1e-4);
  TestAnalyticalStiffArrhenius(rosenbrock_3stage, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_4stage, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_4stage_da, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_6stage_da, 1e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_csc_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_csc_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_in_place_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_in_place_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_in_place_8, 2e-5);
}

TEST(AnalyticalExamples, Branched)
{
  TestAnalyticalBranched(rosenbrock_2stage, 1e-10);
  TestAnalyticalBranched(rosenbrock_3stage);
  TestAnalyticalBranched(rosenbrock_4stage);
  TestAnalyticalBranched(rosenbrock_4stage_da);
  TestAnalyticalBranched(rosenbrock_6stage_da);
  TestAnalyticalBranched(rosenbrock_vector_4);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalBranched(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalBranched(rosenbrock_vector_mozart_in_place_4);
  TestAnalyticalBranched(rosenbrock_vector_mozart_in_place_8);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalStiffBranched(rosenbrock_2stage, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_3stage, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_4stage, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_4stage_da, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_6stage_da, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_csc_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_csc_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_in_place_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_in_place_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_in_place_8, 2e-3);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  TestAnalyticalSurfaceRxn(rosenbrock_2stage, 1e-2);
  TestAnalyticalSurfaceRxn(rosenbrock_3stage, 1e-5);
  TestAnalyticalSurfaceRxn(rosenbrock_4stage, 1e-6);
  TestAnalyticalSurfaceRxn(rosenbrock_4stage_da, 1e-5);
  TestAnalyticalSurfaceRxn(rosenbrock_6stage_da, 1e-7);
  TestAnalyticalSurfaceRxn(rosenbrock_vector_mozart_in_place_8, 1e-5);
}

TEST(AnalyticalExamples, Robertson)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalRobertson(rosenbrock_2stage, 1e-6);
  TestAnalyticalRobertson(rosenbrock_3stage, 1e-6);
  TestAnalyticalRobertson(rosenbrock_4stage, 1e-6);
  TestAnalyticalRobertson(rosenbrock_4stage_da, 1e-6);
  TestAnalyticalRobertson(rosenbrock_6stage_da, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_csc_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_csc_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_in_place_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_in_place_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_in_place_8, 1e-6);
}

TEST(AnalyticalExamples, E5)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalE5(rosenbrock_2stage, 1e-10);
  TestAnalyticalE5(rosenbrock_3stage, 1e-10);
  TestAnalyticalE5(rosenbrock_4stage, 1e-10);
  TestAnalyticalE5(rosenbrock_4stage_da, 1e-10);
  TestAnalyticalE5(rosenbrock_6stage_da, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_csc_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_csc_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_in_place_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_in_place_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_in_place_8, 1e-10);
}

TEST(AnalyticalExamples, Oregonator)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  micm::Real rel_tol = 1e-6;
  TestAnalyticalOregonator(rosenbrock_2stage, rel_tol);
  TestAnalyticalOregonator(rosenbrock_3stage, rel_tol);
  TestAnalyticalOregonator(rosenbrock_4stage, rel_tol);
  TestAnalyticalOregonator(rosenbrock_4stage_da, rel_tol);
  TestAnalyticalOregonator(rosenbrock_6stage_da, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_csc_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_csc_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_in_place_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_in_place_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_in_place_8, rel_tol);
}

TEST(AnalyticalExamples, HIRES)
{
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    GTEST_SKIP() << "Stiff analytical problem is not solvable to the required accuracy in single precision.";
  }

  TestAnalyticalHires(rosenbrock_2stage, 1e-6);
  TestAnalyticalHires(rosenbrock_3stage, 1e-7);
  TestAnalyticalHires(rosenbrock_4stage, 1e-7);
  TestAnalyticalHires(rosenbrock_4stage_da, 1e-6);
  TestAnalyticalHires(rosenbrock_6stage_da, 1e-6);
  TestAnalyticalHires(rosenbrock_vector_doolittle_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_csc_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_csc_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_in_place_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_in_place_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_in_place_8, 1e-7);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
