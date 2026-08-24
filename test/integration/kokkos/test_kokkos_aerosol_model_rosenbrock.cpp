// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "../aerosol_model_policy.hpp"

#include <micm/util/types.hpp>
#include <micm/Kokkos.hpp>

#include <gtest/gtest.h>

#define JUST_ONE_SOLVER

using BuilderType = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>;
using StateType = micm::State<BuilderType::DenseMatrixPolicyType, BuilderType::SparseMatrixPolicyType>;

template<micm::Index L>
using VectorRosenbrock = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::KokkosDenseMatrix<micm::Real, L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

template<micm::Index L>
using VectorRosenbrockDoolittle = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::KokkosDenseMatrix<micm::Real, L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle<micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>>;

template<micm::Index L>
using VectorRosenbrockMozart = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::KokkosDenseMatrix<micm::Real, L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart<micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>>;

template<micm::Index L>
using VectorRosenbrockDolittleCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::KokkosDenseMatrix<micm::Real, L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle<
        micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>>>;

template<micm::Index L>
using VectorRosenbrockMozartCSC = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::KokkosDenseMatrix<micm::Real, L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart<micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>>>;

#ifdef JUST_ONE_SOLVER
auto rosenbrock = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
#else
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

auto rosenbrock_vector_1 = VectorRosenbrock<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_2 = VectorRosenbrock<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_3 = VectorRosenbrock<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_4 = VectorRosenbrock<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

auto rosenbrock_vector_doolittle_1 =
    VectorRosenbrockDoolittle<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_2 =
    VectorRosenbrockDoolittle<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_3 =
    VectorRosenbrockDoolittle<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
auto rosenbrock_vector_doolittle_4 =
    VectorRosenbrockDoolittle<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

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
#endif

TEST(AerosolModelIntegration, StateIncludesStubAerosolModel)
{
#ifdef JUST_ONE_SOLVER
  TestStateIncludesStubAerosolModel(rosenbrock);
#else
  TestStateIncludesStubAerosolModel(rosenbrock_2stage);
  TestStateIncludesStubAerosolModel(rosenbrock_3stage);
  TestStateIncludesStubAerosolModel(rosenbrock_4stage);
  TestStateIncludesStubAerosolModel(rosenbrock_4stage_da);
  TestStateIncludesStubAerosolModel(rosenbrock_6stage_da);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_1);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_2);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_3);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_4);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestStateIncludesStubAerosolModel(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanUpdateStateWithStubAerosolModel)
{
#ifdef JUST_ONE_SOLVER
  TestUpdateStateWithStubAerosolModel(rosenbrock);
#else
  TestUpdateStateWithStubAerosolModel(rosenbrock_2stage);
  TestUpdateStateWithStubAerosolModel(rosenbrock_3stage);
  TestUpdateStateWithStubAerosolModel(rosenbrock_4stage);
  TestUpdateStateWithStubAerosolModel(rosenbrock_4stage_da);
  TestUpdateStateWithStubAerosolModel(rosenbrock_6stage_da);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_1);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_2);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_3);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_4);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestUpdateStateWithStubAerosolModel(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanUpdateMultiCellStateWithStubAerosolModel)
{
#ifdef JUST_ONE_SOLVER
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock);
#else
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_2stage);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_3stage);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_4stage);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_4stage_da);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_6stage_da);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_1);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_2);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_3);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_4);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestUpdateMultiCellStateWithStubAerosolModel(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
#ifdef JUST_ONE_SOLVER
  TestSingleCellForcingWithStubAerosolModel(rosenbrock);
#else
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_2stage);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_3stage);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_4stage);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_4stage_da);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_6stage_da);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_1);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_2);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_3);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_4);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestSingleCellForcingWithStubAerosolModel(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
#ifdef JUST_ONE_SOLVER
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock);
#else
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_2stage);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_3stage);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_4stage);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_4stage_da);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_6stage_da);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_1);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_2);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_3);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_4);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_1);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_2);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_3);
  TestSingleCellJacobianWithStubAerosolModel(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithStubAerosolModel1)
{
#ifdef JUST_ONE_SOLVER
  TestSolveWithStubAerosolModel1(rosenbrock);
#else
  TestSolveWithStubAerosolModel1(rosenbrock_2stage);
  TestSolveWithStubAerosolModel1(rosenbrock_3stage);
  TestSolveWithStubAerosolModel1(rosenbrock_4stage);
  TestSolveWithStubAerosolModel1(rosenbrock_4stage_da);
  TestSolveWithStubAerosolModel1(rosenbrock_6stage_da);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_1);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_2);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_3);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_4);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_1);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_2);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_3);
  TestSolveWithStubAerosolModel1(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithTwoStubAerosolModels)
{
#ifdef JUST_ONE_SOLVER
  TestSolveWithTwoStubAerosolModels(rosenbrock);
#else
  TestSolveWithTwoStubAerosolModels(rosenbrock_2stage, 5e-4);
  TestSolveWithTwoStubAerosolModels(rosenbrock_3stage);
  TestSolveWithTwoStubAerosolModels(rosenbrock_4stage);
  TestSolveWithTwoStubAerosolModels(rosenbrock_4stage_da);
  TestSolveWithTwoStubAerosolModels(rosenbrock_6stage_da);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_1);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_2);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_3);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_4);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_1);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_2);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_3);
  TestSolveWithTwoStubAerosolModels(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithStubAerosolModel1)
{
#ifdef JUST_ONE_SOLVER
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock);
#else
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_2stage);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_3stage);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_4stage);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_4stage_da);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_6stage_da);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_1);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_2);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_3);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_4);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_1);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_2);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_3);
  TestSolveWithStubAerosolModel1MultiCell(rosenbrock_vector_doolittle_4);
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
#endif
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithTwoStubAerosolModels)
{
#ifdef JUST_ONE_SOLVER
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock);
#else
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_2stage, 5e-4);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_3stage);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_4stage);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_4stage_da);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_6stage_da);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_1);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_2);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_3);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_4);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_1);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_2);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_3);
  TestSolveWithTwoStubAerosolModelsMultiCell(rosenbrock_vector_doolittle_4);
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
#endif
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
