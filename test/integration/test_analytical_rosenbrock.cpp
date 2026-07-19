#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/util/types.hpp>

#include <gtest/gtest.h>

using BuilderType = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>;
using StateType = micm::State<BuilderType::DenseMatrixPolicyType, BuilderType::SparseMatrixPolicyType>;

template<micm::Index L>
using VectorRosenbrock = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

using StandardRosenbrockDoolittle = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::Matrix<micm::Real>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionDoolittle>;
using StandardStateTypeDoolittle = micm::State<
    micm::Matrix<micm::Real>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionDoolittle>;

template<micm::Index L>
using VectorRosenbrockDoolittle = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

template<micm::Index L>
using VectorStateTypeDoolittle = micm::State<
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

using StandardRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::Matrix<micm::Real>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionMozart>;

using StandardStateTypeMozart = micm::
    State<micm::Matrix<micm::Real>, micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>, micm::LuDecompositionMozart>;

template<micm::Index L>
using VectorRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<micm::Index L>
using VectorStateTypeMozart = micm::State<
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<micm::Index L>
using VectorStateType =
    micm::State<micm::VectorMatrix<micm::Real, L>, micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

template<micm::Index L>
using VectorRosenbrockDolittleCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle>;

template<micm::Index L>
using VectorStateTypeDoolittleCSC = typename VectorRosenbrockDolittleCSC<L>::StatePolicyType;

template<micm::Index L>
using VectorRosenbrockMozartCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart>;

template<micm::Index L>
using VectorStateTypeMozartCSC = typename VectorRosenbrockMozartCSC<L>::StatePolicyType;

template<micm::Index L>
using VectorRosenbrockDoolittleInPlace = micm::CpuSolverBuilderInPlace<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittleInPlace>;

template<micm::Index L>
using VectorStateTypeDoolittleInPlace = typename VectorRosenbrockDoolittleInPlace<L>::StatePolicyType;

template<micm::Index L>
using VectorRosenbrockMozartInPlace = micm::CpuSolverBuilderInPlace<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace>;

template<micm::Index L>
using VectorStateTypeMozartInPlace = typename VectorRosenbrockMozartInPlace<L>::StatePolicyType;

auto rosenbrock_2stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
auto rosenbrock_3stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_4stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
auto rosenbrock_4stage_da = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
auto rosenbrock_6stage_da = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());

auto rosenbrock_vector_1 = VectorRosenbrock<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_2 = VectorRosenbrock<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_3 = VectorRosenbrock<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_4 = VectorRosenbrock<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_standard_doolittle =
    StandardRosenbrockDoolittle(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_1 =
    VectorRosenbrockDoolittle<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_2 =
    VectorRosenbrockDoolittle<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_3 =
    VectorRosenbrockDoolittle<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_4 =
    VectorRosenbrockDoolittle<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_standard_mozart =
    StandardRosenbrockMozart(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_1 =
    VectorRosenbrockMozart<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_2 =
    VectorRosenbrockMozart<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_3 =
    VectorRosenbrockMozart<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_4 =
    VectorRosenbrockMozart<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_doolittle_csc_1 =
    VectorRosenbrockDolittleCSC<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_csc_2 =
    VectorRosenbrockDolittleCSC<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_csc_3 =
    VectorRosenbrockDolittleCSC<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_csc_4 =
    VectorRosenbrockDolittleCSC<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_mozart_csc_1 =
    VectorRosenbrockMozartCSC<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_csc_2 =
    VectorRosenbrockMozartCSC<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_csc_3 =
    VectorRosenbrockMozartCSC<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_csc_4 =
    VectorRosenbrockMozartCSC<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_doolittle_in_place_1 =
    VectorRosenbrockDoolittleInPlace<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_in_place_2 =
    VectorRosenbrockDoolittleInPlace<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_in_place_3 =
    VectorRosenbrockDoolittleInPlace<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_in_place_4 =
    VectorRosenbrockDoolittleInPlace<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_mozart_in_place_1 =
    VectorRosenbrockMozartInPlace<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_in_place_2 =
    VectorRosenbrockMozartInPlace<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_in_place_3 =
    VectorRosenbrockMozartInPlace<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_in_place_4 =
    VectorRosenbrockMozartInPlace<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

TEST(AnalyticalExamples, Troe)
{
  TestAnalyticalTroe(rosenbrock_2stage);
  TestAnalyticalTroe(rosenbrock_3stage);
  TestAnalyticalTroe(rosenbrock_4stage);
  TestAnalyticalTroe(rosenbrock_4stage_da);
  TestAnalyticalTroe(rosenbrock_6stage_da);
  TestAnalyticalTroe(rosenbrock_vector_1);
  TestAnalyticalTroe(rosenbrock_vector_2);
  TestAnalyticalTroe(rosenbrock_vector_3);
  TestAnalyticalTroe(rosenbrock_vector_4);
  TestAnalyticalTroe(rosenbrock_standard_doolittle);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_1);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_2);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_3);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_4);
  TestAnalyticalTroe(rosenbrock_standard_mozart);
  TestAnalyticalTroe(rosenbrock_vector_mozart_1);
  TestAnalyticalTroe(rosenbrock_vector_mozart_2);
  TestAnalyticalTroe(rosenbrock_vector_mozart_3);
  TestAnalyticalTroe(rosenbrock_vector_mozart_4);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalTroe(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalTroe(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalTroe(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalTroe(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalTroe(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalTroe(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalTroe(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalTroe(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalTroe(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  TestAnalyticalStiffTroe(rosenbrock_2stage);
  TestAnalyticalStiffTroe(rosenbrock_3stage);
  TestAnalyticalStiffTroe(rosenbrock_4stage);
  TestAnalyticalStiffTroe(rosenbrock_4stage_da);
  TestAnalyticalStiffTroe(rosenbrock_6stage_da);
  TestAnalyticalStiffTroe(rosenbrock_vector_1);
  TestAnalyticalStiffTroe(rosenbrock_vector_2);
  TestAnalyticalStiffTroe(rosenbrock_vector_3);
  TestAnalyticalStiffTroe(rosenbrock_vector_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalStiffTroe(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalStiffTroe(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, Photolysis)
{
  TestAnalyticalPhotolysis(rosenbrock_2stage);
  TestAnalyticalPhotolysis(rosenbrock_3stage);
  TestAnalyticalPhotolysis(rosenbrock_4stage);
  TestAnalyticalPhotolysis(rosenbrock_4stage_da);
  TestAnalyticalPhotolysis(rosenbrock_6stage_da);
  TestAnalyticalPhotolysis(rosenbrock_vector_1);
  TestAnalyticalPhotolysis(rosenbrock_vector_2);
  TestAnalyticalPhotolysis(rosenbrock_vector_3);
  TestAnalyticalPhotolysis(rosenbrock_vector_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalPhotolysis(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalPhotolysis(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  TestAnalyticalStiffPhotolysis(rosenbrock_2stage);
  TestAnalyticalStiffPhotolysis(rosenbrock_3stage);
  TestAnalyticalStiffPhotolysis(rosenbrock_4stage);
  TestAnalyticalStiffPhotolysis(rosenbrock_4stage_da);
  TestAnalyticalStiffPhotolysis(rosenbrock_6stage_da);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_1);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_2);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_3);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalStiffPhotolysis(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  TestAnalyticalTernaryChemicalActivation(rosenbrock_2stage);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_3stage);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_4stage);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_4stage_da);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_6stage_da);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_1);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_2);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_3);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_2stage, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_3stage, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_4stage, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_4stage_da, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_6stage_da, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_1, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_2, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_3, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_1, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_2, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_3, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_csc_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_csc_1, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_csc_2, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_csc_3, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_csc_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_1, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_2, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_3, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_doolittle_in_place_4, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_1, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_2, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_3, 2e-3);
  TestAnalyticalStiffTernaryChemicalActivation(rosenbrock_vector_mozart_in_place_4, 2e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  TestAnalyticalTunneling(rosenbrock_2stage, 2e-5);
  TestAnalyticalTunneling(rosenbrock_3stage);
  TestAnalyticalTunneling(rosenbrock_4stage);
  TestAnalyticalTunneling(rosenbrock_4stage_da);
  TestAnalyticalTunneling(rosenbrock_6stage_da);
  TestAnalyticalTunneling(rosenbrock_vector_1);
  TestAnalyticalTunneling(rosenbrock_vector_2);
  TestAnalyticalTunneling(rosenbrock_vector_3);
  TestAnalyticalTunneling(rosenbrock_vector_4);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalTunneling(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalTunneling(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  TestAnalyticalStiffTunneling(rosenbrock_2stage, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_3stage, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_4stage, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_4stage_da, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_6stage_da, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_1, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_2, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_3, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_csc_1, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_csc_2, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_csc_3, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_csc_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_csc_1, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_csc_2, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_csc_3, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_csc_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_in_place_1, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_in_place_2, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_in_place_3, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_doolittle_in_place_4, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_in_place_1, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_in_place_2, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_in_place_3, 1e-4);
  TestAnalyticalStiffTunneling(rosenbrock_vector_mozart_in_place_4, 1e-4);
}

TEST(AnalyticalExamples, Arrhenius)
{
  TestAnalyticalArrhenius(rosenbrock_2stage, 4e-6);
  TestAnalyticalArrhenius(rosenbrock_3stage);
  TestAnalyticalArrhenius(rosenbrock_4stage);
  TestAnalyticalArrhenius(rosenbrock_4stage_da);
  TestAnalyticalArrhenius(rosenbrock_6stage_da);
  TestAnalyticalArrhenius(rosenbrock_vector_1);
  TestAnalyticalArrhenius(rosenbrock_vector_2);
  TestAnalyticalArrhenius(rosenbrock_vector_3);
  TestAnalyticalArrhenius(rosenbrock_vector_4);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalArrhenius(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalArrhenius(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  TestAnalyticalStiffArrhenius(rosenbrock_2stage, 1e-4);
  TestAnalyticalStiffArrhenius(rosenbrock_3stage, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_4stage, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_4stage_da, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_6stage_da, 1e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_1, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_2, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_3, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_csc_1, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_csc_2, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_csc_3, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_csc_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_csc_1, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_csc_2, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_csc_3, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_csc_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_in_place_1, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_in_place_2, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_in_place_3, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_doolittle_in_place_4, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_in_place_1, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_in_place_2, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_in_place_3, 2e-5);
  TestAnalyticalStiffArrhenius(rosenbrock_vector_mozart_in_place_4, 2e-5);
}

TEST(AnalyticalExamples, Branched)
{
  TestAnalyticalBranched(rosenbrock_2stage, 1e-10);
  TestAnalyticalBranched(rosenbrock_3stage);
  TestAnalyticalBranched(rosenbrock_4stage);
  TestAnalyticalBranched(rosenbrock_4stage_da);
  TestAnalyticalBranched(rosenbrock_6stage_da);
  TestAnalyticalBranched(rosenbrock_vector_1);
  TestAnalyticalBranched(rosenbrock_vector_2);
  TestAnalyticalBranched(rosenbrock_vector_3);
  TestAnalyticalBranched(rosenbrock_vector_4);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_csc_1);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_csc_2);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_csc_3);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_csc_4);
  TestAnalyticalBranched(rosenbrock_vector_mozart_csc_1);
  TestAnalyticalBranched(rosenbrock_vector_mozart_csc_2);
  TestAnalyticalBranched(rosenbrock_vector_mozart_csc_3);
  TestAnalyticalBranched(rosenbrock_vector_mozart_csc_4);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_in_place_1);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_in_place_2);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_in_place_3);
  TestAnalyticalBranched(rosenbrock_vector_doolittle_in_place_4);
  TestAnalyticalBranched(rosenbrock_vector_mozart_in_place_1);
  TestAnalyticalBranched(rosenbrock_vector_mozart_in_place_2);
  TestAnalyticalBranched(rosenbrock_vector_mozart_in_place_3);
  TestAnalyticalBranched(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  TestAnalyticalStiffBranched(rosenbrock_2stage, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_3stage, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_4stage, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_4stage_da, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_6stage_da, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_1, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_2, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_3, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_csc_1, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_csc_2, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_csc_3, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_csc_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_csc_1, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_csc_2, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_csc_3, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_csc_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_in_place_1, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_in_place_2, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_in_place_3, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_doolittle_in_place_4, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_in_place_1, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_in_place_2, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_in_place_3, 2e-3);
  TestAnalyticalStiffBranched(rosenbrock_vector_mozart_in_place_4, 2e-3);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  TestAnalyticalSurfaceRxn(rosenbrock_2stage, 1e-2);
  TestAnalyticalSurfaceRxn(rosenbrock_3stage, 1e-5);
  TestAnalyticalSurfaceRxn(rosenbrock_4stage, 1e-6);
  TestAnalyticalSurfaceRxn(rosenbrock_4stage_da, 1e-5);
  TestAnalyticalSurfaceRxn(rosenbrock_6stage_da, 1e-7);
}

TEST(AnalyticalExamples, Robertson)
{
  TestAnalyticalRobertson(rosenbrock_2stage, 1e-6);
  TestAnalyticalRobertson(rosenbrock_3stage, 1e-6);
  TestAnalyticalRobertson(rosenbrock_4stage, 1e-6);
  TestAnalyticalRobertson(rosenbrock_4stage_da, 1e-6);
  TestAnalyticalRobertson(rosenbrock_6stage_da, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_1, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_2, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_3, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_standard_mozart, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_1, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_2, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_3, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_csc_1, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_csc_2, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_csc_3, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_csc_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_csc_1, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_csc_2, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_csc_3, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_csc_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_in_place_1, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_in_place_2, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_in_place_3, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_doolittle_in_place_4, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_in_place_1, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_in_place_2, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_in_place_3, 1e-6);
  TestAnalyticalRobertson(rosenbrock_vector_mozart_in_place_4, 1e-6);
}

TEST(AnalyticalExamples, E5)
{
  TestAnalyticalE5(rosenbrock_2stage, 1e-10);
  TestAnalyticalE5(rosenbrock_3stage, 1e-10);
  TestAnalyticalE5(rosenbrock_4stage, 1e-10);
  TestAnalyticalE5(rosenbrock_4stage_da, 1e-10);
  TestAnalyticalE5(rosenbrock_6stage_da, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_1, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_2, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_3, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_4, 1e-10);
  TestAnalyticalE5(rosenbrock_standard_mozart, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_1, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_2, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_3, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_csc_1, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_csc_2, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_csc_3, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_csc_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_csc_1, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_csc_2, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_csc_3, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_csc_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_in_place_1, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_in_place_2, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_in_place_3, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_doolittle_in_place_4, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_in_place_1, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_in_place_2, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_in_place_3, 1e-10);
  TestAnalyticalE5(rosenbrock_vector_mozart_in_place_4, 1e-10);
}

TEST(AnalyticalExamples, Oregonator)
{
  micm::Real rel_tol = 1e-2;
  TestAnalyticalOregonator(rosenbrock_2stage, rel_tol);
  TestAnalyticalOregonator(rosenbrock_3stage, rel_tol);
  TestAnalyticalOregonator(rosenbrock_4stage, rel_tol);
  TestAnalyticalOregonator(rosenbrock_4stage_da, rel_tol);
  TestAnalyticalOregonator(rosenbrock_6stage_da, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_1, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_2, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_3, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_standard_mozart, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_1, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_2, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_3, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_csc_1, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_csc_2, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_csc_3, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_csc_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_csc_1, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_csc_2, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_csc_3, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_csc_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_in_place_1, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_in_place_2, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_in_place_3, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_doolittle_in_place_4, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_in_place_1, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_in_place_2, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_in_place_3, rel_tol);
  TestAnalyticalOregonator(rosenbrock_vector_mozart_in_place_4, rel_tol);
}

TEST(AnalyticalExamples, HIRES)
{
  TestAnalyticalHires(rosenbrock_2stage, 1e-6);
  TestAnalyticalHires(rosenbrock_3stage, 1e-7);
  TestAnalyticalHires(rosenbrock_4stage, 1e-7);
  TestAnalyticalHires(rosenbrock_4stage_da, 1e-6);
  TestAnalyticalHires(rosenbrock_6stage_da, 1e-6);
  TestAnalyticalHires(rosenbrock_vector_doolittle_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_4, 1e-7);
  TestAnalyticalHires(rosenbrock_standard_mozart, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_csc_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_csc_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_csc_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_csc_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_csc_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_csc_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_csc_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_csc_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_in_place_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_in_place_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_in_place_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_doolittle_in_place_4, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_in_place_1, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_in_place_2, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_in_place_3, 1e-7);
  TestAnalyticalHires(rosenbrock_vector_mozart_in_place_4, 1e-7);
}
