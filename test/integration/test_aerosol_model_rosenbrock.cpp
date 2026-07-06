// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "aerosol_model_policy.hpp"

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

template<std::size_t L>
using VectorRosenbrockDoolittle = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

using StandardRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::Matrix<double>,
    micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorRosenbrockDolittleCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorRosenbrockMozartCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart>;

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

TEST(AerosolModelIntegration, StateIncludesStubAerosolModel)
{
  TestStateIncludesStubAerosolModel(rosenbrock_2stage);
  TestStateIncludesStubAerosolModel(rosenbrock_3stage);
  TestStateIncludesStubAerosolModel(rosenbrock_4stage);
  TestStateIncludesStubAerosolModel(rosenbrock_4stage_da);
  TestStateIncludesStubAerosolModel(rosenbrock_6stage_da);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_1);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_2);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_3);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_4);
  TestStateIncludesStubAerosolModel(rosenbrock_standard_doolittle);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_4);
  TestStateIncludesStubAerosolModel(rosenbrock_standard_mozart);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_1);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_2);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_3);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_4);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_csc_1);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_csc_2);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_csc_3);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_csc_4);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_csc_1);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_csc_2);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_csc_3);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanUpdateStateWithStubAerosolModel)
{
  TestUpdateStateWithStubAerosolModel(rosenbrock_2stage);
  TestUpdateStateWithStubAerosolModel(rosenbrock_3stage);
  TestUpdateStateWithStubAerosolModel(rosenbrock_4stage);
  TestUpdateStateWithStubAerosolModel(rosenbrock_4stage_da);
  TestUpdateStateWithStubAerosolModel(rosenbrock_6stage_da);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_1);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_2);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_3);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_4);
  TestUpdateStateWithStubAerosolModel(rosenbrock_standard_doolittle);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_4);
  TestUpdateStateWithStubAerosolModel(rosenbrock_standard_mozart);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_1);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_2);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_3);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_4);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_1);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_2);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_3);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_4);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_1);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_2);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_3);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanUpdateMultiCellStateWithStubAerosolModel)
{
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_2stage);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_3stage);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_4stage);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_4stage_da);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_6stage_da);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_1);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_2);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_3);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_4);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_standard_doolittle);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_4);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_standard_mozart);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_1);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_2);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_3);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_4);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_1);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_2);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_3);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_csc_4);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_1);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_2);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_3);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_2stage);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_3stage);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_4stage);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_4stage_da);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_6stage_da);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_1);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_2);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_3);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_4);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_standard_doolittle);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_4);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_standard_mozart);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_1);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_2);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_3);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_4);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_csc_1);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_csc_2);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_csc_3);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_csc_4);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_csc_1);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_csc_2);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_csc_3);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_2stage);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_3stage);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_4stage);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_4stage_da);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_6stage_da);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_1);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_2);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_3);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_4);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_standard_doolittle);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_4);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_standard_mozart);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_1);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_2);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_3);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_4);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_csc_1);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_csc_2);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_csc_3);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_csc_4);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_csc_1);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_csc_2);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_csc_3);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithStubAerosolModel1)
{
  TestSolveWithStubAerosolModel1(rosenbrock_2stage);
  TestSolveWithStubAerosolModel1(rosenbrock_3stage);
  TestSolveWithStubAerosolModel1(rosenbrock_4stage);
  TestSolveWithStubAerosolModel1(rosenbrock_4stage_da);
  TestSolveWithStubAerosolModel1(rosenbrock_6stage_da);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_1);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_2);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_3);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_4);
  TestSolveWithStubAerosolModel1(rosenbrock_standard_doolittle);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_1);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_2);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_3);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_4);
  TestSolveWithStubAerosolModel1(rosenbrock_standard_mozart);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_1);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_2);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_3);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_4);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_csc_1);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_csc_2);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_csc_3);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_csc_4);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_csc_1);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_csc_2);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_csc_3);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithTwoStubAerosolModels)
{
  TestSolveWithTwoStubAerosolModels(rosenbrock_2stage, 5e-4);
  TestSolveWithTwoStubAerosolModels(rosenbrock_3stage);
  TestSolveWithTwoStubAerosolModels(rosenbrock_4stage);
  TestSolveWithTwoStubAerosolModels(rosenbrock_4stage_da);
  TestSolveWithTwoStubAerosolModels(rosenbrock_6stage_da);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_1);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_2);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_3);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_4);
  TestSolveWithTwoStubAerosolModels(rosenbrock_standard_doolittle);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_1);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_2);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_3);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_4);
  TestSolveWithTwoStubAerosolModels(rosenbrock_standard_mozart);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_1);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_2);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_3);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_4);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_csc_1);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_csc_2);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_csc_3);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_csc_4);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_csc_1);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_csc_2);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_csc_3);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithStubAerosolModel1)
{
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_2stage);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_3stage);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_4stage);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_4stage_da);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_6stage_da);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_1);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_2);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_3);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_4);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_standard_doolittle);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_1);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_2);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_3);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_4);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_standard_mozart);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_1);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_2);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_3);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_4);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_csc_1);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_csc_2);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_csc_3);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_csc_4);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_csc_1);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_csc_2);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_csc_3);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithTwoStubAerosolModels)
{
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_2stage, 5e-4);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_3stage);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_4stage);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_4stage_da);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_6stage_da);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_1);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_2);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_3);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_4);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_standard_doolittle);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_1);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_2);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_3);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_4);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_standard_mozart);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_1);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_2);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_3);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_4);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_csc_1);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_csc_2);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_csc_3);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_csc_4);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_csc_1);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_csc_2);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_csc_3);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_mozart_csc_4);
}
