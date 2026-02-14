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
  test_state_includes_stub_aerosol_model(rosenbrock_2stage);
  test_state_includes_stub_aerosol_model(rosenbrock_3stage);
  test_state_includes_stub_aerosol_model(rosenbrock_4stage);
  test_state_includes_stub_aerosol_model(rosenbrock_4stage_da);
  test_state_includes_stub_aerosol_model(rosenbrock_6stage_da);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_1);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_2);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_3);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_4);
  test_state_includes_stub_aerosol_model(rosenbrock_standard_doolittle);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_1);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_2);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_3);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_4);
  test_state_includes_stub_aerosol_model(rosenbrock_standard_mozart);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_1);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_2);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_3);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_4);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_csc_1);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_csc_2);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_csc_3);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_doolittle_csc_4);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_csc_1);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_csc_2);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_csc_3);
  test_state_includes_stub_aerosol_model(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanUpdateStateWithStubAerosolModel)
{
  test_update_state_with_stub_aerosol_model(rosenbrock_2stage);
  test_update_state_with_stub_aerosol_model(rosenbrock_3stage);
  test_update_state_with_stub_aerosol_model(rosenbrock_4stage);
  test_update_state_with_stub_aerosol_model(rosenbrock_4stage_da);
  test_update_state_with_stub_aerosol_model(rosenbrock_6stage_da);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_1);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_2);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_3);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_4);
  test_update_state_with_stub_aerosol_model(rosenbrock_standard_doolittle);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_1);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_2);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_3);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_4);
  test_update_state_with_stub_aerosol_model(rosenbrock_standard_mozart);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_1);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_2);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_3);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_4);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_1);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_2);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_3);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_4);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_1);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_2);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_3);
  test_update_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanUpdateMultiCellStateWithStubAerosolModel)
{
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_2stage);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_3stage);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_4stage);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_4stage_da);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_6stage_da);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_1);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_2);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_3);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_4);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_standard_doolittle);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_1);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_2);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_3);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_4);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_standard_mozart);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_1);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_2);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_3);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_4);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_1);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_2);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_3);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_4);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_1);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_2);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_3);
  test_update_multi_cell_state_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_2stage);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_3stage);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_4stage);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_4stage_da);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_6stage_da);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_1);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_2);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_3);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_4);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_standard_doolittle);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_1);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_2);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_3);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_4);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_standard_mozart);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_1);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_2);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_3);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_4);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_1);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_2);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_3);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_4);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_1);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_2);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_3);
  test_single_cell_forcing_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_2stage);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_3stage);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_4stage);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_4stage_da);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_6stage_da);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_1);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_2);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_3);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_4);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_standard_doolittle);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_1);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_2);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_3);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_4);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_standard_mozart);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_1);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_2);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_3);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_4);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_1);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_2);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_3);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_doolittle_csc_4);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_1);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_2);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_3);
  test_single_cell_jacobian_with_stub_aerosol_model(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithStubAerosolModel1)
{
  test_solve_with_stub_aerosol_model_1(rosenbrock_2stage);
  test_solve_with_stub_aerosol_model_1(rosenbrock_3stage);
  test_solve_with_stub_aerosol_model_1(rosenbrock_4stage);
  test_solve_with_stub_aerosol_model_1(rosenbrock_4stage_da);
  test_solve_with_stub_aerosol_model_1(rosenbrock_6stage_da);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_1);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_2);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_3);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_4);
  test_solve_with_stub_aerosol_model_1(rosenbrock_standard_doolittle);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_1);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_2);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_3);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_4);
  test_solve_with_stub_aerosol_model_1(rosenbrock_standard_mozart);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_1);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_2);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_3);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_4);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_csc_1);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_csc_2);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_csc_3);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_doolittle_csc_4);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_csc_1);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_csc_2);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_csc_3);
  test_solve_with_stub_aerosol_model_1(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveSingleGridCellWithTwoStubAerosolModels)
{
  test_solve_with_two_stub_aerosol_models(rosenbrock_2stage, 5e-4);
  test_solve_with_two_stub_aerosol_models(rosenbrock_3stage);
  test_solve_with_two_stub_aerosol_models(rosenbrock_4stage);
  test_solve_with_two_stub_aerosol_models(rosenbrock_4stage_da);
  test_solve_with_two_stub_aerosol_models(rosenbrock_6stage_da);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_1);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_2);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_3);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_4);
  test_solve_with_two_stub_aerosol_models(rosenbrock_standard_doolittle);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_1);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_2);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_3);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_4);
  test_solve_with_two_stub_aerosol_models(rosenbrock_standard_mozart);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_1);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_2);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_3);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_4);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_csc_1);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_csc_2);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_csc_3);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_doolittle_csc_4);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_csc_1);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_csc_1);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_csc_2);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_csc_3);
  test_solve_with_two_stub_aerosol_models(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithStubAerosolModel1)
{
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_2stage);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_3stage);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_4stage);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_4stage_da);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_6stage_da);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_1);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_2);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_3);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_4);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_standard_doolittle);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_1);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_2);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_3);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_4);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_standard_mozart);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_1);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_2);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_3);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_4);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_csc_1);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_csc_2);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_csc_3);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_doolittle_csc_4);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_csc_1);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_csc_2);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_csc_3);
  test_solve_with_stub_aerosol_model_1_multi_cell(rosenbrock_vector_mozart_csc_4);
}

TEST(AerosolModelIntegration, CanSolveMultiGridCellWithTwoStubAerosolModels)
{
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_2stage, 5e-4);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_3stage);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_4stage);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_4stage_da);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_6stage_da);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_1);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_2);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_3);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_4);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_standard_doolittle);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_1);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_2);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_3);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_4);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_standard_mozart);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_1);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_2);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_3);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_4);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_csc_1);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_csc_2);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_csc_3);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_doolittle_csc_4);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_csc_1);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_csc_1);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_csc_2);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_csc_3);
  test_solve_with_two_stub_aerosol_models_multi_cell(rosenbrock_vector_mozart_csc_4);
}

