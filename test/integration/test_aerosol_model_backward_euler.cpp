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
  test_state_includes_stub_aerosol_model(backward_euler);
  test_state_includes_stub_aerosol_model(backward_euler_vector_1);
  test_state_includes_stub_aerosol_model(backward_euler_vector_2);
  test_state_includes_stub_aerosol_model(backward_euler_vector_3);
  test_state_includes_stub_aerosol_model(backward_euler_vector_4);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_1);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_2);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_3);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_4);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_1);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_2);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_3);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_4);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_csc_1);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_csc_2);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_csc_3);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_csc_4);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_csc_1);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_csc_2);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_csc_3);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_csc_4);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_in_place_1);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_in_place_2);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_in_place_3);
  test_state_includes_stub_aerosol_model(backward_euler_vector_doolittle_in_place_4);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_in_place_1);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_in_place_2);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_in_place_3);
  test_state_includes_stub_aerosol_model(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanUpdateStateWithStubAerosolModel)
{
  test_update_state_with_stub_aerosol_model(backward_euler);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_1);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_2);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_3);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_4);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_1);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_2);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_3);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_4);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_1);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_2);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_3);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_4);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_1);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_2);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_3);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_4);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_1);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_2);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_3);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_4);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_1);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_2);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_3);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_4);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_1);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_2);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_3);
  test_update_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanUpdateMultiCellStateWithStubAerosolModel)
{
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_1);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_2);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_3);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_4);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_1);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_2);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_3);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_4);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_1);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_2);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_3);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_4);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_1);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_2);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_3);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_4);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_1);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_2);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_3);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_csc_4);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_1);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_2);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_3);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_4);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_1);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_2);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_3);
  test_update_multi_cell_state_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_1);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_2);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_3);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_4);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_1);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_2);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_3);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_4);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_1);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_2);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_3);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_4);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_1);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_2);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_3);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_4);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_csc_1);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_csc_2);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_csc_3);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_csc_4);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_1);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_2);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_3);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_4);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_1);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_2);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_3);
  test_single_cell_forcing_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_1);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_2);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_3);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_4);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_1);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_2);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_3);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_4);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_1);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_2);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_3);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_4);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_1);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_2);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_3);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_csc_4);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_csc_1);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_csc_2);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_csc_3);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_csc_4);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_1);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_2);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_3);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_doolittle_in_place_4);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_1);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_2);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_3);
  test_single_cell_jacobian_with_stub_aerosol_model(backward_euler_vector_mozart_in_place_4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithStubAerosolModel1)
{
  test_solve_with_stub_aerosol_model_1(backward_euler, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_1, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_2, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_3, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_4, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_1, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_2, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_3, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_4, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_1, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_2, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_3, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_4, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_csc_1, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_csc_2, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_csc_3, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_csc_4, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_csc_1, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_csc_2, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_csc_3, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_csc_4, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_in_place_1, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_in_place_2, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_in_place_3, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_doolittle_in_place_4, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_in_place_1, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_in_place_2, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_in_place_3, 2e-4);
  test_solve_with_stub_aerosol_model_1(backward_euler_vector_mozart_in_place_4, 2e-4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithTwoStubAerosolModels)
{
  test_solve_with_two_stub_aerosol_models(backward_euler, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_1, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_2, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_3, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_4, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_1, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_2, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_3, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_4, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_1, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_2, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_3, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_4, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_csc_1, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_csc_2, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_csc_3, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_csc_4, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_csc_1, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_csc_2, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_csc_3, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_csc_4, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_in_place_1, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_in_place_2, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_in_place_3, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_doolittle_in_place_4, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_in_place_1, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_in_place_2, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_in_place_3, 1e-1);
  test_solve_with_two_stub_aerosol_models(backward_euler_vector_mozart_in_place_4, 1e-1);
}
