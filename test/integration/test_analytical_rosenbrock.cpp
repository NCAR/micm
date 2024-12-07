#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

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
using StandardStateTypeDoolittle = micm::State<micm::Matrix<double>, micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>, micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorRosenbrockDoolittle = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorStateTypeDoolittle =
    micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>, micm::LuDecompositionDoolittle>;

using StandardRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::Matrix<double>,
    micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionMozart>;

using StandardStateTypeMozart = micm::State<micm::Matrix<double>, micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>, micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateTypeMozart =
    micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>, micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorStateType =
    micm::State<micm::VectorMatrix<double, L>, micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

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

auto rosenbrock_standard_doolittle = StandardRosenbrockDoolittle(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_1 = VectorRosenbrockDoolittle<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_2 = VectorRosenbrockDoolittle<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_3 = VectorRosenbrockDoolittle<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_4 = VectorRosenbrockDoolittle<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_standard_mozart = StandardRosenbrockMozart(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_1 = VectorRosenbrockMozart<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_2 = VectorRosenbrockMozart<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_3 = VectorRosenbrockMozart<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_mozart_4 = VectorRosenbrockMozart<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe(rosenbrock_2stage);
  test_analytical_troe(rosenbrock_3stage);
  test_analytical_troe(rosenbrock_4stage);
  test_analytical_troe(rosenbrock_4stage_da);
  test_analytical_troe(rosenbrock_6stage_da);
  test_analytical_troe<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_troe<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_troe<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_troe<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
  test_analytical_troe<StandardRosenbrockDoolittle, StandardStateTypeDoolittle>(rosenbrock_standard_doolittle);
  test_analytical_troe<VectorRosenbrockDoolittle<1>, VectorStateTypeDoolittle<1>>(rosenbrock_vector_doolittle_1);
  test_analytical_troe<VectorRosenbrockDoolittle<2>, VectorStateTypeDoolittle<2>>(rosenbrock_vector_doolittle_2);
  test_analytical_troe<VectorRosenbrockDoolittle<3>, VectorStateTypeDoolittle<3>>(rosenbrock_vector_doolittle_3);
  test_analytical_troe<VectorRosenbrockDoolittle<4>, VectorStateTypeDoolittle<4>>(rosenbrock_vector_doolittle_4);
  test_analytical_troe<StandardRosenbrockMozart, StandardStateTypeMozart>(rosenbrock_standard_mozart);
  test_analytical_troe<VectorRosenbrockMozart<1>, VectorStateTypeMozart<1>>(rosenbrock_vector_mozart_1);
  test_analytical_troe<VectorRosenbrockMozart<2>, VectorStateTypeMozart<2>>(rosenbrock_vector_mozart_2);
  test_analytical_troe<VectorRosenbrockMozart<3>, VectorStateTypeMozart<3>>(rosenbrock_vector_mozart_3);
  test_analytical_troe<VectorRosenbrockMozart<4>, VectorStateTypeMozart<4>>(rosenbrock_vector_mozart_4);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe(rosenbrock_2stage);
  test_analytical_stiff_troe(rosenbrock_3stage);
  test_analytical_stiff_troe(rosenbrock_4stage);
  test_analytical_stiff_troe(rosenbrock_4stage_da);
  test_analytical_stiff_troe(rosenbrock_6stage_da);
  test_analytical_stiff_troe<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_stiff_troe<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_stiff_troe<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_stiff_troe<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis(rosenbrock_2stage);
  test_analytical_photolysis(rosenbrock_3stage);
  test_analytical_photolysis(rosenbrock_4stage);
  test_analytical_photolysis(rosenbrock_4stage_da);
  test_analytical_photolysis(rosenbrock_6stage_da);
  test_analytical_photolysis<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_photolysis<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_photolysis<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_photolysis<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis(rosenbrock_2stage);
  test_analytical_stiff_photolysis(rosenbrock_3stage);
  test_analytical_stiff_photolysis(rosenbrock_4stage);
  test_analytical_stiff_photolysis(rosenbrock_4stage_da);
  test_analytical_stiff_photolysis(rosenbrock_6stage_da);
  test_analytical_stiff_photolysis<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_stiff_photolysis<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_stiff_photolysis<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_stiff_photolysis<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation(rosenbrock_2stage);
  test_analytical_ternary_chemical_activation(rosenbrock_3stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage);
  test_analytical_ternary_chemical_activation(rosenbrock_4stage_da);
  test_analytical_ternary_chemical_activation(rosenbrock_6stage_da);
  test_analytical_ternary_chemical_activation<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_ternary_chemical_activation<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_ternary_chemical_activation<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_ternary_chemical_activation<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_2stage, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_3stage, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_4stage_da, 2e-3);
  test_analytical_stiff_ternary_chemical_activation(rosenbrock_6stage_da, 2e-3);
  test_analytical_stiff_ternary_chemical_activation<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1, 2e-3);
  test_analytical_stiff_ternary_chemical_activation<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2, 2e-3);
  test_analytical_stiff_ternary_chemical_activation<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3, 2e-3);
  test_analytical_stiff_ternary_chemical_activation<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4, 2e-3);
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling(rosenbrock_2stage, 2e-5);
  test_analytical_tunneling(rosenbrock_3stage);
  test_analytical_tunneling(rosenbrock_4stage);
  test_analytical_tunneling(rosenbrock_4stage_da);
  test_analytical_tunneling(rosenbrock_6stage_da);
  test_analytical_tunneling<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_tunneling<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_tunneling<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_tunneling<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_3stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_4stage_da, 1e-4);
  test_analytical_stiff_tunneling(rosenbrock_6stage_da, 1e-4);
  test_analytical_stiff_tunneling<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1, 1e-4);
  test_analytical_stiff_tunneling<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2, 1e-4);
  test_analytical_stiff_tunneling<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3, 1e-4);
  test_analytical_stiff_tunneling<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4, 1e-4);
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius(rosenbrock_2stage, 4e-6);
  test_analytical_arrhenius(rosenbrock_3stage);
  test_analytical_arrhenius(rosenbrock_4stage);
  test_analytical_arrhenius(rosenbrock_4stage_da);
  test_analytical_arrhenius(rosenbrock_6stage_da);
  test_analytical_arrhenius<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_arrhenius<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_arrhenius<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_arrhenius<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius(rosenbrock_2stage, 1e-4);
  test_analytical_stiff_arrhenius(rosenbrock_3stage, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_4stage_da, 2e-5);
  test_analytical_stiff_arrhenius(rosenbrock_6stage_da, 1e-5);
  test_analytical_stiff_arrhenius<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1, 2e-5);
  test_analytical_stiff_arrhenius<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2, 2e-5);
  test_analytical_stiff_arrhenius<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3, 2e-5);
  test_analytical_stiff_arrhenius<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4, 2e-5);
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched(rosenbrock_2stage, 1e-10);
  test_analytical_branched(rosenbrock_3stage);
  test_analytical_branched(rosenbrock_4stage);
  test_analytical_branched(rosenbrock_4stage_da);
  test_analytical_branched(rosenbrock_6stage_da);
  test_analytical_branched<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1);
  test_analytical_branched<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2);
  test_analytical_branched<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3);
  test_analytical_branched<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched(rosenbrock_2stage, 2e-3);
  test_analytical_stiff_branched(rosenbrock_3stage, 2e-3);
  test_analytical_stiff_branched(rosenbrock_4stage, 2e-3);
  test_analytical_stiff_branched(rosenbrock_4stage_da, 2e-3);
  test_analytical_stiff_branched(rosenbrock_6stage_da, 2e-3);
  test_analytical_stiff_branched<VectorRosenbrock<1>, VectorStateType<1>>(rosenbrock_vector_1, 2e-3);
  test_analytical_stiff_branched<VectorRosenbrock<2>, VectorStateType<2>>(rosenbrock_vector_2, 2e-3);
  test_analytical_stiff_branched<VectorRosenbrock<3>, VectorStateType<3>>(rosenbrock_vector_3, 2e-3);
  test_analytical_stiff_branched<VectorRosenbrock<4>, VectorStateType<4>>(rosenbrock_vector_4, 2e-3);
}

TEST(AnalyticalExamples, Robertson)
{
  auto rosenbrock_solver = [](auto params)
  {
    params.relative_tolerance_ = 1e-10;
    params.absolute_tolerance_ = std::vector<double>(3, params.relative_tolerance_ * 1e-2);
    return BuilderType(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_robertson<BuilderType, StateType>(solver, 1e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_robertson<BuilderType, StateType>(solver, 1e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_robertson<BuilderType, StateType>(solver, 1e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_robertson<BuilderType, StateType>(solver, 1e-1);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_robertson<BuilderType, StateType>(solver, 1e-1);

  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  params.relative_tolerance_ = 1e-10;
  params.absolute_tolerance_ = std::vector<double>(3, params.relative_tolerance_ * 1e-2);

  auto standard_dootlittle = StandardRosenbrockDoolittle(params);
  test_analytical_robertson<StandardRosenbrockDoolittle, StandardStateTypeDoolittle>(standard_dootlittle, 1e-1);
  auto vector_1_dootlittle = VectorRosenbrockDoolittle<1>(params);
  test_analytical_robertson<VectorRosenbrockDoolittle<1>, VectorStateTypeDoolittle<1>>(vector_1_dootlittle, 1e-1);
  auto vector_2_dootlittle = VectorRosenbrockDoolittle<2>(params);
  test_analytical_robertson<VectorRosenbrockDoolittle<2>, VectorStateTypeDoolittle<2>>(vector_2_dootlittle, 1e-1);
  auto vector_3_dootlittle = VectorRosenbrockDoolittle<3>(params);
  test_analytical_robertson<VectorRosenbrockDoolittle<3>, VectorStateTypeDoolittle<3>>(vector_3_dootlittle, 1e-1);
  auto vector_4_dootlittle = VectorRosenbrockDoolittle<4>(params);
  test_analytical_robertson<VectorRosenbrockDoolittle<4>, VectorStateTypeDoolittle<4>>(vector_4_dootlittle, 1e-1);
  auto standard_mozart = StandardRosenbrockMozart(params);
  test_analytical_robertson<StandardRosenbrockMozart, StandardStateTypeMozart>(standard_mozart, 1e-1);
  auto vector_1_mozart = VectorRosenbrockMozart<1>(params);
  test_analytical_robertson<VectorRosenbrockMozart<1>, VectorStateTypeMozart<1>>(vector_1_mozart, 1e-1);
  auto vector_2_mozart = VectorRosenbrockMozart<2>(params);
  test_analytical_robertson<VectorRosenbrockMozart<2>, VectorStateTypeMozart<2>>(vector_2_mozart, 1e-1);
  auto vector_3_mozart = VectorRosenbrockMozart<3>(params);
  test_analytical_robertson<VectorRosenbrockMozart<3>, VectorStateTypeMozart<3>>(vector_3_mozart, 1e-1);
  auto vector_4_mozart = VectorRosenbrockMozart<4>(params);
  test_analytical_robertson<VectorRosenbrockMozart<4>, VectorStateTypeMozart<4>>(vector_4_mozart, 1e-1);
}

TEST(AnalyticalExamples, E5)
{
  auto rosenbrock_solver = [](auto params)
  {
    params.relative_tolerance_ = 1e-13;
    params.absolute_tolerance_ = std::vector<double>(6, 1e-17);
    // this paper https://archimede.uniba.it/~testset/report/e5.pdf
    // says that the first variable should have a much looser tolerance than the other species
    params.absolute_tolerance_[0] = 1e-7;
    // these last two aren't actually provided values and we don't care how they behave
    params.absolute_tolerance_[4] = 1e-7;
    params.absolute_tolerance_[5] = 1e-7;
    return BuilderType(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_e5<BuilderType, StateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_e5<BuilderType, StateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_e5<BuilderType, StateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_e5<BuilderType, StateType>(solver, 1e-3);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_e5<BuilderType, StateType>(solver, 1e-3);

  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  params.relative_tolerance_ = 1e-13;
  params.absolute_tolerance_ = std::vector<double>(6, 1e-17);
  params.absolute_tolerance_[0] = 1e-7;
  params.absolute_tolerance_[4] = 1e-7;
  params.absolute_tolerance_[5] = 1e-7;

  auto standard_dootlittle = StandardRosenbrockDoolittle(params);
  test_analytical_e5<StandardRosenbrockDoolittle, StandardStateTypeDoolittle>(standard_dootlittle, 1e-3);
  auto vector_1_dootlittle = VectorRosenbrockDoolittle<1>(params);
  test_analytical_e5<VectorRosenbrockDoolittle<1>, VectorStateTypeDoolittle<1>>(vector_1_dootlittle, 1e-3);
  auto vector_2_dootlittle = VectorRosenbrockDoolittle<2>(params);
  test_analytical_e5<VectorRosenbrockDoolittle<2>, VectorStateTypeDoolittle<2>>(vector_2_dootlittle, 1e-3);
  auto vector_3_dootlittle = VectorRosenbrockDoolittle<3>(params);
  test_analytical_e5<VectorRosenbrockDoolittle<3>, VectorStateTypeDoolittle<3>>(vector_3_dootlittle, 1e-3);
  auto vector_4_dootlittle = VectorRosenbrockDoolittle<4>(params);
  test_analytical_e5<VectorRosenbrockDoolittle<4>, VectorStateTypeDoolittle<4>>(vector_4_dootlittle, 1e-3);
  auto standard_mozart = StandardRosenbrockMozart(params);
  test_analytical_e5<StandardRosenbrockMozart, StandardStateTypeMozart>(standard_mozart, 1e-3);
  auto vector_1_mozart = VectorRosenbrockMozart<1>(params);
  test_analytical_e5<VectorRosenbrockMozart<1>, VectorStateTypeMozart<1>>(vector_1_mozart, 1e-3);
  auto vector_2_mozart = VectorRosenbrockMozart<2>(params);
  test_analytical_e5<VectorRosenbrockMozart<2>, VectorStateTypeMozart<2>>(vector_2_mozart, 1e-3);
  auto vector_3_mozart = VectorRosenbrockMozart<3>(params);
  test_analytical_e5<VectorRosenbrockMozart<3>, VectorStateTypeMozart<3>>(vector_3_mozart, 1e-3);
  auto vector_4_mozart = VectorRosenbrockMozart<4>(params);
  test_analytical_e5<VectorRosenbrockMozart<4>, VectorStateTypeMozart<4>>(vector_4_mozart, 1e-3);

}

TEST(AnalyticalExamples, Oregonator)
{
  auto rosenbrock_solver = [](auto params)
  {
    // anything below 1e-6 is too strict for the Oregonator
    params.relative_tolerance_ = 1e-8;
    params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-6);
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_oregonator(solver, 4e-6);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_oregonator(solver, 4e-6);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_oregonator(solver, 4e-6);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_oregonator(solver, 4e-6);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_oregonator(solver, 4e-6);

  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  params.relative_tolerance_ = 1e-8;
  params.absolute_tolerance_ = std::vector<double>(5, params.relative_tolerance_ * 1e-6);

  auto standard_dootlittle = StandardRosenbrockDoolittle(params);
  test_analytical_oregonator<StandardRosenbrockDoolittle, StandardStateTypeDoolittle>(standard_dootlittle, 4e-6);
  auto vector_1_dootlittle = VectorRosenbrockDoolittle<1>(params);
  test_analytical_oregonator<VectorRosenbrockDoolittle<1>, VectorStateTypeDoolittle<1>>(vector_1_dootlittle, 4e-6);
  auto vector_2_dootlittle = VectorRosenbrockDoolittle<2>(params);
  test_analytical_oregonator<VectorRosenbrockDoolittle<2>, VectorStateTypeDoolittle<2>>(vector_2_dootlittle, 4e-6);
  auto vector_3_dootlittle = VectorRosenbrockDoolittle<3>(params);
  test_analytical_oregonator<VectorRosenbrockDoolittle<3>, VectorStateTypeDoolittle<3>>(vector_3_dootlittle, 4e-6);
  auto vector_4_dootlittle = VectorRosenbrockDoolittle<4>(params);
  test_analytical_oregonator<VectorRosenbrockDoolittle<4>, VectorStateTypeDoolittle<4>>(vector_4_dootlittle, 4e-6);
  auto standard_mozart = StandardRosenbrockMozart(params);
  test_analytical_oregonator<StandardRosenbrockMozart, StandardStateTypeMozart>(standard_mozart, 4e-6);
  auto vector_1_mozart = VectorRosenbrockMozart<1>(params);
  test_analytical_oregonator<VectorRosenbrockMozart<1>, VectorStateTypeMozart<1>>(vector_1_mozart, 4e-6);
  auto vector_2_mozart = VectorRosenbrockMozart<2>(params);
  test_analytical_oregonator<VectorRosenbrockMozart<2>, VectorStateTypeMozart<2>>(vector_2_mozart, 4e-6);
  auto vector_3_mozart = VectorRosenbrockMozart<3>(params);
  test_analytical_oregonator<VectorRosenbrockMozart<3>, VectorStateTypeMozart<3>>(vector_3_mozart, 4e-6);
  auto vector_4_mozart = VectorRosenbrockMozart<4>(params);
  test_analytical_oregonator<VectorRosenbrockMozart<4>, VectorStateTypeMozart<4>>(vector_4_mozart, 4e-6);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn(rosenbrock_2stage, 1e-2);
  test_analytical_surface_rxn(rosenbrock_3stage, 1e-5);
  test_analytical_surface_rxn(rosenbrock_4stage, 1e-6);
  test_analytical_surface_rxn(rosenbrock_4stage_da, 1e-5);
  test_analytical_surface_rxn(rosenbrock_6stage_da, 1e-7);
}

TEST(AnalyticalExamples, HIRES)
{
  auto rosenbrock_solver = [](auto params)
  {
    params.relative_tolerance_ = 1e-6;
    params.absolute_tolerance_ = std::vector<double>(8, params.relative_tolerance_ * 1e-2);
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params);
  };

  auto solver = rosenbrock_solver(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  test_analytical_hires(solver, 1e-6);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_hires(solver, 1e-7);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  test_analytical_hires(solver, 1e-7);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_hires(solver, 1e-6);

  solver = rosenbrock_solver(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  test_analytical_hires(solver, 1e-6);

  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  params.relative_tolerance_ = 1e-6;
  params.absolute_tolerance_ = std::vector<double>(8, params.relative_tolerance_ * 1e-2);
  
  auto standard_dootlittle = StandardRosenbrockDoolittle(params);
  test_analytical_hires<StandardRosenbrockDoolittle, StandardStateTypeDoolittle>(standard_dootlittle, 1e-7);
  auto vector_1_dootlittle = VectorRosenbrockDoolittle<1>(params);
  test_analytical_hires<VectorRosenbrockDoolittle<1>, VectorStateTypeDoolittle<1>>(vector_1_dootlittle, 1e-7);
  auto vector_2_dootlittle = VectorRosenbrockDoolittle<2>(params);
  test_analytical_hires<VectorRosenbrockDoolittle<2>, VectorStateTypeDoolittle<2>>(vector_2_dootlittle, 1e-7);
  auto vector_3_dootlittle = VectorRosenbrockDoolittle<3>(params);
  test_analytical_hires<VectorRosenbrockDoolittle<3>, VectorStateTypeDoolittle<3>>(vector_3_dootlittle, 1e-7);
  auto vector_4_dootlittle = VectorRosenbrockDoolittle<4>(params);
  test_analytical_hires<VectorRosenbrockDoolittle<4>, VectorStateTypeDoolittle<4>>(vector_4_dootlittle, 1e-7);
  auto standard_mozart = StandardRosenbrockMozart(params);
  test_analytical_hires<StandardRosenbrockMozart, StandardStateTypeMozart>(standard_mozart, 1e-7);
  auto vector_1_mozart = VectorRosenbrockMozart<1>(params);
  test_analytical_hires<VectorRosenbrockMozart<1>, VectorStateTypeMozart<1>>(vector_1_mozart, 1e-7);
  auto vector_2_mozart = VectorRosenbrockMozart<2>(params);
  test_analytical_hires<VectorRosenbrockMozart<2>, VectorStateTypeMozart<2>>(vector_2_mozart, 1e-7);
  auto vector_3_mozart = VectorRosenbrockMozart<3>(params);
  test_analytical_hires<VectorRosenbrockMozart<3>, VectorStateTypeMozart<3>>(vector_3_mozart, 1e-7);
  auto vector_4_mozart = VectorRosenbrockMozart<4>(params);
  test_analytical_hires<VectorRosenbrockMozart<4>, VectorStateTypeMozart<4>>(vector_4_mozart, 1e-7);
}
