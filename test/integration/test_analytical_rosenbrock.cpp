#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <gtest/gtest.h>

using BuilderType = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>;
using StateType = micm::State<BuilderType::DenseMatrixPolicyType, BuilderType::SparseMatrixPolicyType>;

template<std::size_t L>
using VectorRosenbrock = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

using StandardRosenbrockDoolittle = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::Matrix<double>,
    micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionDoolittle>;
using StandardStateTypeDoolittle = micm::State<
    micm::Matrix<double>,
    micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorRosenbrockDoolittle = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorStateTypeDoolittle = micm::State<
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

using StandardRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::Matrix<double>,
    micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionMozart>;

using StandardStateTypeMozart = micm::
    State<micm::Matrix<double>, micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>, micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateTypeMozart = micm::State<
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateType =
    micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

template<std::size_t L>
using VectorRosenbrockDolittleCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorStateTypeDoolittleCSC = typename VectorRosenbrockDolittleCSC<L>::StatePolicyType;

template<std::size_t L>
using VectorRosenbrockMozartCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateTypeMozartCSC = typename VectorRosenbrockMozartCSC<L>::StatePolicyType;

template<std::size_t L>
using VectorRosenbrockDoolittleInPlace = micm::CpuSolverBuilderInPlace<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittleInPlace>;

template<std::size_t L>
using VectorStateTypeDoolittleInPlace = typename VectorRosenbrockDoolittleInPlace<L>::StatePolicyType;

template<std::size_t L>
using VectorRosenbrockMozartInPlace = micm::CpuSolverBuilderInPlace<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace>;

template<std::size_t L>
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
  test_analytical_troe(rosenbrock_2stage);
  test_analytical_troe(rosenbrock_3stage);
  test_analytical_troe(rosenbrock_4stage);
  test_analytical_troe(rosenbrock_4stage_da);
  test_analytical_troe(rosenbrock_6stage_da);
  test_analytical_troe(rosenbrock_vector_1);
  test_analytical_troe(rosenbrock_vector_2);
  test_analytical_troe(rosenbrock_vector_3);
  test_analytical_troe(rosenbrock_vector_4);
  test_analytical_troe(rosenbrock_standard_doolittle);
  test_analytical_troe(rosenbrock_vector_doolittle_1);
  test_analytical_troe(rosenbrock_vector_doolittle_2);
  test_analytical_troe(rosenbrock_vector_doolittle_3);
  test_analytical_troe(rosenbrock_vector_doolittle_4);
  test_analytical_troe(rosenbrock_standard_mozart);
  test_analytical_troe(rosenbrock_vector_mozart_1);
  test_analytical_troe(rosenbrock_vector_mozart_2);
  test_analytical_troe(rosenbrock_vector_mozart_3);
  test_analytical_troe(rosenbrock_vector_mozart_4);
  test_analytical_troe(rosenbrock_vector_doolittle_csc_1);
  test_analytical_troe(rosenbrock_vector_doolittle_csc_2);
  test_analytical_troe(rosenbrock_vector_doolittle_csc_3);
  test_analytical_troe(rosenbrock_vector_doolittle_csc_4);
  test_analytical_troe(rosenbrock_vector_mozart_csc_1);
  test_analytical_troe(rosenbrock_vector_mozart_csc_2);
  test_analytical_troe(rosenbrock_vector_mozart_csc_3);
  test_analytical_troe(rosenbrock_vector_mozart_csc_4);
  test_analytical_troe(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_troe(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_troe(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_troe(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_troe(rosenbrock_vector_mozart_in_place_1);
  test_analytical_troe(rosenbrock_vector_mozart_in_place_2);
  test_analytical_troe(rosenbrock_vector_mozart_in_place_3);
  test_analytical_troe(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(rosenbrock_2stage);
  test_analytical_stiff_troe(rosenbrock_3stage);
  test_analytical_stiff_troe(rosenbrock_4stage);
  test_analytical_stiff_troe(rosenbrock_4stage_da);
  test_analytical_stiff_troe(rosenbrock_6stage_da);
  test_analytical_stiff_troe(rosenbrock_vector_1);
  test_analytical_stiff_troe(rosenbrock_vector_2);
  test_analytical_stiff_troe(rosenbrock_vector_3);
  test_analytical_stiff_troe(rosenbrock_vector_4);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_csc_1);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_csc_2);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_csc_3);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_csc_4);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_csc_1);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_csc_2);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_csc_3);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_csc_4);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_stiff_troe(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_in_place_1);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_in_place_2);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_in_place_3);
  test_analytical_stiff_troe(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(rosenbrock_2stage);
  test_analytical_photolysis(rosenbrock_3stage);
  test_analytical_photolysis(rosenbrock_4stage);
  test_analytical_photolysis(rosenbrock_4stage_da);
  test_analytical_photolysis(rosenbrock_6stage_da);
  test_analytical_photolysis(rosenbrock_vector_1);
  test_analytical_photolysis(rosenbrock_vector_2);
  test_analytical_photolysis(rosenbrock_vector_3);
  test_analytical_photolysis(rosenbrock_vector_4);
  test_analytical_photolysis(rosenbrock_vector_doolittle_csc_1);
  test_analytical_photolysis(rosenbrock_vector_doolittle_csc_2);
  test_analytical_photolysis(rosenbrock_vector_doolittle_csc_3);
  test_analytical_photolysis(rosenbrock_vector_doolittle_csc_4);
  test_analytical_photolysis(rosenbrock_vector_mozart_csc_1);
  test_analytical_photolysis(rosenbrock_vector_mozart_csc_2);
  test_analytical_photolysis(rosenbrock_vector_mozart_csc_3);
  test_analytical_photolysis(rosenbrock_vector_mozart_csc_4);
  test_analytical_photolysis(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_photolysis(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_photolysis(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_photolysis(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_photolysis(rosenbrock_vector_mozart_in_place_1);
  test_analytical_photolysis(rosenbrock_vector_mozart_in_place_2);
  test_analytical_photolysis(rosenbrock_vector_mozart_in_place_3);
  test_analytical_photolysis(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(rosenbrock_2stage);
  test_analytical_stiff_photolysis(rosenbrock_3stage);
  test_analytical_stiff_photolysis(rosenbrock_4stage);
  test_analytical_stiff_photolysis(rosenbrock_4stage_da);
  test_analytical_stiff_photolysis(rosenbrock_6stage_da);
  test_analytical_stiff_photolysis(rosenbrock_vector_1);
  test_analytical_stiff_photolysis(rosenbrock_vector_2);
  test_analytical_stiff_photolysis(rosenbrock_vector_3);
  test_analytical_stiff_photolysis(rosenbrock_vector_4);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_csc_1);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_csc_2);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_csc_3);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_csc_4);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_csc_1);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_csc_2);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_csc_3);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_csc_4);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_stiff_photolysis(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_in_place_1);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_in_place_2);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_in_place_3);
  test_analytical_stiff_photolysis(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(rosenbrock_2stage);
  test_analytical_ternary_chemical_activation(rosenbrock_3stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage_da);
  test_analytical_ternary_chemical_activation(rosenbrock_6stage_da);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_1);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_2);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_3);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_4);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_1);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_2);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_3);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_4);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_csc_1);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_csc_2);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_csc_3);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_csc_4);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_1);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_2);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_3);
  test_analytical_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_2stage, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_3stage, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage_da, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_6stage_da, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_1, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_2, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_3, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_4, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_1, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_2, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_3, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_csc_4, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_csc_1, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_csc_2, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_csc_3, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_csc_4, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_1, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_2, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_3, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_doolittle_in_place_4, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_1, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_2, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_3, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_vector_mozart_in_place_4, 2e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(rosenbrock_2stage, 2e-5);
  test_analytical_tunneling(rosenbrock_3stage);
  test_analytical_tunneling(rosenbrock_4stage);
  test_analytical_tunneling(rosenbrock_4stage_da);
  test_analytical_tunneling(rosenbrock_6stage_da);
  test_analytical_tunneling(rosenbrock_vector_1);
  test_analytical_tunneling(rosenbrock_vector_2);
  test_analytical_tunneling(rosenbrock_vector_3);
  test_analytical_tunneling(rosenbrock_vector_4);
  test_analytical_tunneling(rosenbrock_vector_doolittle_csc_1);
  test_analytical_tunneling(rosenbrock_vector_doolittle_csc_2);
  test_analytical_tunneling(rosenbrock_vector_doolittle_csc_3);
  test_analytical_tunneling(rosenbrock_vector_doolittle_csc_4);
  test_analytical_tunneling(rosenbrock_vector_mozart_csc_1);
  test_analytical_tunneling(rosenbrock_vector_mozart_csc_2);
  test_analytical_tunneling(rosenbrock_vector_mozart_csc_3);
  test_analytical_tunneling(rosenbrock_vector_mozart_csc_4);
  test_analytical_tunneling(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_tunneling(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_tunneling(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_tunneling(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_tunneling(rosenbrock_vector_mozart_in_place_1);
  test_analytical_tunneling(rosenbrock_vector_mozart_in_place_2);
  test_analytical_tunneling(rosenbrock_vector_mozart_in_place_3);
  test_analytical_tunneling(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_6stage_da, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_1, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_2, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_3, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_4, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_csc_1, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_csc_2, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_csc_3, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_csc_4, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_csc_1, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_csc_2, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_csc_3, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_csc_4, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_in_place_1, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_in_place_2, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_in_place_3, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_doolittle_in_place_4, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_in_place_1, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_in_place_2, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_in_place_3, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_vector_mozart_in_place_4, 1e-4);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(rosenbrock_2stage, 4e-6);
  test_analytical_arrhenius(rosenbrock_3stage);
  test_analytical_arrhenius(rosenbrock_4stage);
  test_analytical_arrhenius(rosenbrock_4stage_da);
  test_analytical_arrhenius(rosenbrock_6stage_da);
  test_analytical_arrhenius(rosenbrock_vector_1);
  test_analytical_arrhenius(rosenbrock_vector_2);
  test_analytical_arrhenius(rosenbrock_vector_3);
  test_analytical_arrhenius(rosenbrock_vector_4);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_csc_1);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_csc_2);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_csc_3);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_csc_4);
  test_analytical_arrhenius(rosenbrock_vector_mozart_csc_1);
  test_analytical_arrhenius(rosenbrock_vector_mozart_csc_2);
  test_analytical_arrhenius(rosenbrock_vector_mozart_csc_3);
  test_analytical_arrhenius(rosenbrock_vector_mozart_csc_4);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_arrhenius(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_arrhenius(rosenbrock_vector_mozart_in_place_1);
  test_analytical_arrhenius(rosenbrock_vector_mozart_in_place_2);
  test_analytical_arrhenius(rosenbrock_vector_mozart_in_place_3);
  test_analytical_arrhenius(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_1, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_2, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_3, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_4, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_csc_1, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_csc_2, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_csc_3, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_csc_4, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_csc_1, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_csc_2, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_csc_3, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_csc_4, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_in_place_1, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_in_place_2, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_in_place_3, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_doolittle_in_place_4, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_in_place_1, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_in_place_2, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_in_place_3, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_vector_mozart_in_place_4, 2e-5);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched(rosenbrock_2stage, 1e-10);
  test_analytical_branched(rosenbrock_3stage);
  test_analytical_branched(rosenbrock_4stage);
  test_analytical_branched(rosenbrock_4stage_da);
  test_analytical_branched(rosenbrock_6stage_da);
  test_analytical_branched(rosenbrock_vector_1);
  test_analytical_branched(rosenbrock_vector_2);
  test_analytical_branched(rosenbrock_vector_3);
  test_analytical_branched(rosenbrock_vector_4);
  test_analytical_branched(rosenbrock_vector_doolittle_csc_1);
  test_analytical_branched(rosenbrock_vector_doolittle_csc_2);
  test_analytical_branched(rosenbrock_vector_doolittle_csc_3);
  test_analytical_branched(rosenbrock_vector_doolittle_csc_4);
  test_analytical_branched(rosenbrock_vector_mozart_csc_1);
  test_analytical_branched(rosenbrock_vector_mozart_csc_2);
  test_analytical_branched(rosenbrock_vector_mozart_csc_3);
  test_analytical_branched(rosenbrock_vector_mozart_csc_4);
  test_analytical_branched(rosenbrock_vector_doolittle_in_place_1);
  test_analytical_branched(rosenbrock_vector_doolittle_in_place_2);
  test_analytical_branched(rosenbrock_vector_doolittle_in_place_3);
  test_analytical_branched(rosenbrock_vector_doolittle_in_place_4);
  test_analytical_branched(rosenbrock_vector_mozart_in_place_1);
  test_analytical_branched(rosenbrock_vector_mozart_in_place_2);
  test_analytical_branched(rosenbrock_vector_mozart_in_place_3);
  test_analytical_branched(rosenbrock_vector_mozart_in_place_4);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(rosenbrock_2stage, 2e-3);
  test_analytical_stiff_branched(rosenbrock_3stage, 2e-3);
  test_analytical_stiff_branched(rosenbrock_4stage, 2e-3);
  test_analytical_stiff_branched(rosenbrock_4stage_da, 2e-3);
  test_analytical_stiff_branched(rosenbrock_6stage_da, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_1, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_2, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_3, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_4, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_csc_1, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_csc_2, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_csc_3, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_csc_4, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_csc_1, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_csc_2, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_csc_3, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_csc_4, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_in_place_1, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_in_place_2, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_in_place_3, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_doolittle_in_place_4, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_in_place_1, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_in_place_2, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_in_place_3, 2e-3);
  test_analytical_stiff_branched(rosenbrock_vector_mozart_in_place_4, 2e-3);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(rosenbrock_2stage, 1e-2);
  test_analytical_surface_rxn(rosenbrock_3stage, 1e-5);
  test_analytical_surface_rxn(rosenbrock_4stage, 1e-6);
  test_analytical_surface_rxn(rosenbrock_4stage_da, 1e-5);
  test_analytical_surface_rxn(rosenbrock_6stage_da, 1e-7);
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson(rosenbrock_2stage, 1e-1);
  test_analytical_robertson(rosenbrock_3stage, 1e-1);
  test_analytical_robertson(rosenbrock_4stage, 1e-1);
  test_analytical_robertson(rosenbrock_4stage_da, 1e-1);
  test_analytical_robertson(rosenbrock_6stage_da, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_1, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_2, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_3, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_4, 1e-1);
  test_analytical_robertson(rosenbrock_standard_mozart, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_1, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_2, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_3, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_4, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_csc_1, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_csc_2, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_csc_3, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_csc_4, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_csc_1, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_csc_2, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_csc_3, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_csc_4, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_in_place_1, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_in_place_2, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_in_place_3, 1e-1);
  test_analytical_robertson(rosenbrock_vector_doolittle_in_place_4, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_in_place_1, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_in_place_2, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_in_place_3, 1e-1);
  test_analytical_robertson(rosenbrock_vector_mozart_in_place_4, 1e-1);
}

TEST(AnalyticalExamples, E5)
{
  test_analytical_e5(rosenbrock_2stage, 1e-3);
  test_analytical_e5(rosenbrock_3stage, 1e-3);
  test_analytical_e5(rosenbrock_4stage, 1e-3);
  test_analytical_e5(rosenbrock_4stage_da, 1e-3);
  test_analytical_e5(rosenbrock_6stage_da, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_1, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_2, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_3, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_4, 1e-3);
  test_analytical_e5(rosenbrock_standard_mozart, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_1, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_2, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_3, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_4, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_csc_1, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_csc_2, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_csc_3, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_csc_4, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_csc_1, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_csc_2, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_csc_3, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_csc_4, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_in_place_1, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_in_place_2, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_in_place_3, 1e-3);
  test_analytical_e5(rosenbrock_vector_doolittle_in_place_4, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_in_place_1, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_in_place_2, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_in_place_3, 1e-3);
  test_analytical_e5(rosenbrock_vector_mozart_in_place_4, 1e-3);
}

TEST(AnalyticalExamples, Oregonator)
{
  test_analytical_oregonator(rosenbrock_2stage, 4e-6);
  test_analytical_oregonator(rosenbrock_3stage, 4e-6);
  test_analytical_oregonator(rosenbrock_4stage, 4e-6);
  test_analytical_oregonator(rosenbrock_4stage_da, 4.5e-6);
  test_analytical_oregonator(rosenbrock_6stage_da, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_1, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_2, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_3, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_4, 4e-6);
  test_analytical_oregonator(rosenbrock_standard_mozart, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_1, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_2, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_3, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_4, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_csc_1, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_csc_2, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_csc_3, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_csc_4, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_csc_1, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_csc_2, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_csc_3, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_csc_4, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_in_place_1, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_in_place_2, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_in_place_3, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_doolittle_in_place_4, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_in_place_1, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_in_place_2, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_in_place_3, 4e-6);
  test_analytical_oregonator(rosenbrock_vector_mozart_in_place_4, 4e-6);
}

TEST(AnalyticalExamples, HIRES)
{
  test_analytical_hires(rosenbrock_2stage, 1e-6);
  test_analytical_hires(rosenbrock_3stage, 1e-7);
  test_analytical_hires(rosenbrock_4stage, 1e-7);
  test_analytical_hires(rosenbrock_4stage_da, 1e-6);
  test_analytical_hires(rosenbrock_6stage_da, 1e-6);
  test_analytical_hires(rosenbrock_vector_doolittle_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_4, 1e-7);
  test_analytical_hires(rosenbrock_standard_mozart, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_4, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_4, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_4, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_csc_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_csc_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_csc_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_csc_4, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_csc_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_csc_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_csc_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_csc_4, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_in_place_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_in_place_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_in_place_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_doolittle_in_place_4, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_in_place_1, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_in_place_2, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_in_place_3, 1e-7);
  test_analytical_hires(rosenbrock_vector_mozart_in_place_4, 1e-7);
}
