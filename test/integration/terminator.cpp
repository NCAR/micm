#include "terminator.hpp"

#include <gtest/gtest.h>

#include <micm/process/process.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void RunTerminatorTest(std::size_t number_of_grid_cells)
{
  TestTerminator<MatrixPolicy, micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>>(
      [&](const micm::System& s, const std::vector<micm::Process>& p)
          -> micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>
      {
        auto solver_params = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells, true);
        solver_params.absolute_tolerance_ = 1.0e-20;
        solver_params.relative_tolerance_ = 1.0e-8;
        solver_params.max_number_of_steps_ = 100000;
        return micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>{ s, p, solver_params };
      },
      number_of_grid_cells);
  TestTerminator<MatrixPolicy, micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>>(
      [&](const micm::System& s, const std::vector<micm::Process>& p)
          -> micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>
      {
        auto solver_params = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells, true);
        solver_params.absolute_tolerance_ = 1.0e-20;
        solver_params.relative_tolerance_ = 1.0e-8;
        solver_params.max_number_of_steps_ = 100000;
        solver_params.check_singularity_ = true;
        return micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>{ s, p, solver_params };
      },
      number_of_grid_cells);
}

TEST(RosenbrockSolver, Terminator)
{
  RunTerminatorTest<micm::Matrix, SparseMatrixTest, micm::LinearSolver<double, SparseMatrixTest>>(2);
  RunTerminatorTest<micm::Matrix, SparseMatrixTest, micm::LinearSolver<double, SparseMatrixTest>>(2);
  RunTerminatorTest<micm::Matrix, SparseMatrixTest, micm::LinearSolver<double, SparseMatrixTest>>(3);
  RunTerminatorTest<micm::Matrix, SparseMatrixTest, micm::LinearSolver<double, SparseMatrixTest>>(4);
}

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

TEST(RosenbrockSolver, VectorTerminator)
{
  RunTerminatorTest<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(1);
  RunTerminatorTest<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(4);
  RunTerminatorTest<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(3);
  RunTerminatorTest<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(2);
}