// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "../aerosol_model_policy.hpp"

#include <micm/Kokkos.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#define JUST_ONE_SOLVER

using BuilderType = micm::KokkosSolverBuilder<micm::BackwardEulerSolverParameters>;
using StateType = micm::State<BuilderType::DenseMatrixPolicyType, BuilderType::SparseMatrixPolicyType>;

template<micm::Index L>
using VectorBackwardEuler = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::KokkosDenseMatrix<micm::Real, L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

template<micm::Index L>
using VectorBackwardEulerMozartInPlace = micm::CpuSolverBuilder<
    micm::BackwardEulerSolverParameters,
    micm::KokkosDenseMatrix<micm::Real, L>,
    micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>,
    micm::LuDecompositionMozartInPlace<micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>>;

auto backward_euler = micm::KokkosSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());

#ifndef JUST_ONE_SOLVER
auto backward_euler_vector_1 = VectorBackwardEuler<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_2 = VectorBackwardEuler<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_3 = VectorBackwardEuler<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_4 = VectorBackwardEuler<4>(micm::BackwardEulerSolverParameters());

auto backward_euler_vector_mozart_in_place_1 = VectorBackwardEulerMozartInPlace<1>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_2 = VectorBackwardEulerMozartInPlace<2>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_3 = VectorBackwardEulerMozartInPlace<3>(micm::BackwardEulerSolverParameters());
auto backward_euler_vector_mozart_in_place_4 = VectorBackwardEulerMozartInPlace<4>(micm::BackwardEulerSolverParameters());
#endif

TEST(AerosolModelIntegration, StateIncludesStubAerosolModel)
{
  TestStateIncludesStubAerosolModel(backward_euler);
#ifndef JUST_ONE_SOLVER
  TestStateIncludesStubAerosolModel(backward_euler_vector_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_4);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestStateIncludesStubAerosolModel(backward_euler_vector_mozart_in_place_4);
#endif
}

TEST(AerosolModelIntegration, CanUpdateStateWithStubAerosolModel)
{
  TestUpdateStateWithStubAerosolModel(backward_euler);
#ifndef JUST_ONE_SOLVER
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_4);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestUpdateStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
#endif
}

TEST(AerosolModelIntegration, CanUpdateMultiCellStateWithStubAerosolModel)
{
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler);
#ifndef JUST_ONE_SOLVER
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_4);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestUpdateMultiCellStateWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
#endif
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
  TestSingleCellForcingWithStubAerosolModel(backward_euler);
#ifndef JUST_ONE_SOLVER
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_4);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestSingleCellForcingWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
#endif
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
  TestSingleCellJacobianWithStubAerosolModel(backward_euler);
#ifndef JUST_ONE_SOLVER
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_4);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_1);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_2);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_3);
  TestSingleCellJacobianWithStubAerosolModel(backward_euler_vector_mozart_in_place_4);
#endif
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithStubAerosolModel1)
{
  TestSolveWithStubAerosolModel1(backward_euler, 2e-4);
#ifndef JUST_ONE_SOLVER
  TestSolveWithStubAerosolModel1(backward_euler_vector_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_4, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_1, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_2, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_3, 2e-4);
  TestSolveWithStubAerosolModel1(backward_euler_vector_mozart_in_place_4, 2e-4);
#endif
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithTwoStubAerosolModels)
{
  TestSolveWithTwoStubAerosolModels(backward_euler, 1e-1);
#ifndef JUST_ONE_SOLVER
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_4, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_1, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_2, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_3, 1e-1);
  TestSolveWithTwoStubAerosolModels(backward_euler_vector_mozart_in_place_4, 1e-1);
#endif
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithStubAerosolModel1)
{
  TestSolveWithStubAerosolModel1MultiCell(backward_euler, 2e-4);
#ifndef JUST_ONE_SOLVER
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_4, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_1, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_2, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_3, 2e-4);
  TestSolveWithStubAerosolModel1MultiCell(backward_euler_vector_mozart_in_place_4, 2e-4);
#endif
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithTwoStubAerosolModels)
{
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler, 1e-1);
#ifndef JUST_ONE_SOLVER
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_4, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_1, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_2, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_3, 1e-1);
  TestSolveWithTwoStubAerosolModelsMultiCell(backward_euler_vector_mozart_in_place_4, 1e-1);
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
