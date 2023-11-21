#pragma once
#include <gtest/gtest.h>

#include <functional>
#include <micm/solver/cuda_linear_solver.hpp>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>

#include "test_linear_solver_policy.hpp"

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;
template<class T>
using Group10000VectorMatrix = micm::VectorMatrix<T, 10000>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;
template<class T>
using Group10000SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<10000>>;

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
std::vector<double> linearSolverGenerator(
    const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver,
    std::size_t number_of_blocks)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy<double>::create(10).number_of_blocks(number_of_blocks).initial_value(1.0e-30);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.with_element(i, j);

  SparseMatrixPolicy<double> A(builder);
  MatrixPolicy<double> b(number_of_blocks, 10, 0.0);
  MatrixPolicy<double> x(number_of_blocks, 10, 100.0);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = get_double();

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      b[i_block][i] = get_double();

  LinearSolverPolicy solver = create_linear_solver(A, 1.0e-30);
  auto lu = micm::LuDecomposition::GetLUMatrices(A, 1.0e-30);
  auto lower_matrix = std::move(lu.first);
  auto upper_matrix = std::move(lu.second);
  solver.Factor(A, lower_matrix, upper_matrix);
  solver.template Solve<MatrixPolicy>(b, x, lower_matrix, upper_matrix);
  return x.AsVector();
}

// bit to bit variation between CPU and GPU result with randomMatrixVectorOrdering
void gpuValidation()
{
  std::vector<double> cpu_x = linearSolverGenerator<
      Group10000VectorMatrix,
      Group10000SparseVectorMatrix,
      micm::LinearSolver<double, Group10000SparseVectorMatrix>>(
      [](const Group10000SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group10000SparseVectorMatrix> {
        return micm::LinearSolver<double, Group10000SparseVectorMatrix>{ matrix, initial_value };
      },
      10000);

  std::vector<double> gpu_x = linearSolverGenerator<
      Group10000VectorMatrix,
      Group10000SparseVectorMatrix,
      micm::CudaLinearSolver<double, Group10000SparseVectorMatrix>>(
      [](const Group10000SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group10000SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group10000SparseVectorMatrix>{ matrix, initial_value };
      },
      10000);

  for (int i = 0; i < cpu_x.size(); i++)
  {
    EXPECT_EQ(cpu_x, gpu_x);
  }
}

TEST(CudaLinearSolver, DenseMatrixVectorOrdering)
{
  testDenseMatrix<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::CudaLinearSolver<double, Group1SparseVectorMatrix, micm::CudaLuDecomposition>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1SparseVectorMatrix, micm::CudaLuDecomposition> {
        return micm::CudaLinearSolver<double, Group1SparseVectorMatrix, micm::CudaLuDecomposition>{ matrix, initial_value };
      });
}

TEST(CudaLinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::CudaLinearSolver<double, Group1SparseVectorMatrix>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group1SparseVectorMatrix>{ matrix, initial_value };
      },
      1);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::CudaLinearSolver<double, Group2SparseVectorMatrix>>(
      [](const Group2SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group2SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group2SparseVectorMatrix>{ matrix, initial_value };
      },
      2);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::CudaLinearSolver<double, Group3SparseVectorMatrix>>(
      [](const Group3SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group3SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group3SparseVectorMatrix>{ matrix, initial_value };
      },
      3);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::CudaLinearSolver<double, Group4SparseVectorMatrix>>(
      [](const Group4SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group4SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group4SparseVectorMatrix>{ matrix, initial_value };
      },
      4);
  gpuValidation();
}

TEST(CudaLinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::CudaLinearSolver<double, Group1SparseVectorMatrix>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group1SparseVectorMatrix>{ matrix, initial_value };
      },
      1);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::CudaLinearSolver<double, Group2SparseVectorMatrix>>(
      [](const Group2SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group2SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group2SparseVectorMatrix>{ matrix, initial_value };
      },
      2);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::CudaLinearSolver<double, Group3SparseVectorMatrix>>(
      [](const Group3SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group3SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group3SparseVectorMatrix>{ matrix, initial_value };
      },
      3);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::CudaLinearSolver<double, Group4SparseVectorMatrix>>(
      [](const Group4SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group4SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group4SparseVectorMatrix>{ matrix, initial_value };
      },
      4);
}
