// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "aerosol_model_policy.hpp"

#include <gtest/gtest.h>

using BuilderType = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>;
using StateType = micm::State<BuilderType::DenseMatrixPolicyType, BuilderType::SparseMatrixPolicyType>;

template<std::size_t L>
using VectorBackwardEuler = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

template<std::size_t L>
using VectorBackwardEulerDoolittle = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorBackwardEulerMozart = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorBackwardEulerDolittleCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionDoolittle>;

template<std::size_t L>
using VectorBackwardEulerMozartCSC = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<L>>,
    micm::LuDecompositionMozart>;

template<std::size_t L>
using VectorBackwardEulerDoolittleInPlace = micm::CpuSolverBuilderInPlace<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionDoolittleInPlace>;

template<std::size_t L>
using VectorBackwardEulerMozartInPlace = micm::CpuSolverBuilderInPlace<
    micm::BackwardEulerSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace>;

auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());

auto backward_euler_vector_1 = VectorBackwardEuler<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_2 = VectorBackwardEuler<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_3 = VectorBackwardEuler<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());

auto backward_euler_vector_doolittle_1 = VectorBackwardEulerDoolittle<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_2 = VectorBackwardEulerDoolittle<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_3 = VectorBackwardEulerDoolittle<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_4 = VectorBackwardEulerDoolittle<4>(micm::BackwardEulerSolverParameters());

auto backward_euler_vector_mozart_1 = VectorBackwardEulerMozart<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_2 = VectorBackwardEulerMozart<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_3 = VectorBackwardEulerMozart<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_4 = VectorBackwardEulerMozart<4>(micm::BackwardEulerSolverParameters());

auto backward_euler_vector_doolittle_csc_1 = VectorBackwardEulerDolittleCSC<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_2 = VectorBackwardEulerDolittleCSC<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_3 = VectorBackwardEulerDolittleCSC<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_doolittle_csc_4 = VectorBackwardEulerDolittleCSC<4>(micm::BackwardEulerSolverParameters());

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

TEST(AerosolModelIntegration, StateIncludesStubAerosolModel)
{
  TestStateIncludesStubAerosolModel(backward_euler);
  TestStateIncludesStubAerosolModel(backward_euler_vector_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_4);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_4);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_4);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_csc_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_csc_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_csc_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_csc_4);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_csc_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_csc_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_csc_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_csc_4);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_in_place_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_in_place_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_in_place_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_doolittle_in_place_4);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanUpdateStateWithStubAerosolModel)
{
  TestUpdateStateWithStubAerosolModel(backward_euler);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_4);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_4);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_4);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_4);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_csc_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_csc_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_csc_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_csc_4);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_4);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanUpdateMultiCellStateWithStubAerosolModel)
{
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_4);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_4);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_4);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_csc_4);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_csc_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_csc_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_csc_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_csc_4);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_doolittle_in_place_4);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
  TestSingleCellForcingWithStubAerosolModel(backward_euler);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_4);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_4);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_4);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_csc_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_csc_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_csc_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_csc_4);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_csc_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_csc_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_csc_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_csc_4);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_in_place_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_in_place_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_in_place_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_doolittle_in_place_4);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
  TestSingleCellJacobianWithStubAerosolModel(backward_euler);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_4);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_4);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_4);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_csc_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_csc_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_csc_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_csc_4);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_csc_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_csc_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_csc_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_csc_4);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_in_place_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_in_place_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_in_place_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_doolittle_in_place_4);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithStubAerosolModel1)
{
  TestSolveWithStubAerosolModel1(backward_euler, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_4, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_4, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_4, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_csc_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_csc_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_csc_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_csc_4, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_csc_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_csc_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_csc_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_csc_4, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_in_place_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_in_place_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_in_place_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_doolittle_in_place_4, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_4, 2e-4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithTwoStubAerosolModels)
{
  TestSolveWithTwoStubAerosolModels(backward_euler, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_4, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_4, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_4, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_csc_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_csc_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_csc_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_csc_4, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_csc_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_csc_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_csc_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_csc_4, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_in_place_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_in_place_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_in_place_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_doolittle_in_place_4, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_4, 1e-1);
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithStubAerosolModel1)
{
  TestSolveWithStubAerosolModel1MultiCell(backward_euler, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_4, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_4, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_4, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_csc_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_csc_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_csc_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_csc_4, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_csc_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_csc_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_csc_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_csc_4, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_in_place_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_in_place_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_in_place_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_doolittle_in_place_4, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_4, 2e-4);
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithTwoStubAerosolModels)
{
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_4, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_4, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_4, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_csc_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_csc_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_csc_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_csc_4, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_csc_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_csc_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_csc_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_csc_4, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_in_place_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_in_place_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_in_place_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_doolittle_in_place_4, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_4, 1e-1);
}
