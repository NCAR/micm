#pragma once
#include "test_linear_solver_policy.hpp"

#include <micm/util/cuda_dense_matrix.hpp>
#include <micm/util/cuda_sparse_matrix.hpp>
#include <micm/solver/cuda_linear_solver.hpp>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <functional>
#include <random>

template<class T>
using Group10000VectorMatrix = micm::VectorMatrix<T, 10000>;

template<class T>
using Group1CudaDenseMatrix = micm::CudaDenseMatrix<T, 1>;
template<class T>
using Group20CudaDenseMatrix = micm::CudaDenseMatrix<T, 20>;
template<class T>
using Group300CudaDenseMatrix = micm::CudaDenseMatrix<T, 300>;
template<class T>
using Group4000CudaDenseMatrix = micm::CudaDenseMatrix<T, 4000>;
template<class T>
using Group10000CudaDenseMatrix = micm::CudaDenseMatrix<T, 10000>;

template<class T>
using Group10000SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<10000>>;

template<class T>
using Group1CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group20CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<20>>;
template<class T>
using Group300CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<300>>;
template<class T>
using Group4000CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<4000>>;
template<class T>
using Group10000CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<10000>>;

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

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(A);
  CopyToDeviceDense<double, MatrixPolicy>(b);
  CopyToDeviceDense<double, MatrixPolicy>(x);

  LinearSolverPolicy solver = create_linear_solver(A, 1.0e-30);
  std::pair<SparseMatrixPolicy<double>, SparseMatrixPolicy<double>> lu = micm::LuDecomposition::GetLUMatrices<double, SparseMatrixPolicy>(A, 1.0e-30);
  SparseMatrixPolicy<double> lower_matrix = std::move(lu.first);
  SparseMatrixPolicy<double> upper_matrix = std::move(lu.second);
  
  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(lower_matrix);
  CopyToDeviceSparse<double, SparseMatrixPolicy>(upper_matrix);

  solver.Factor(A, lower_matrix, upper_matrix);
  solver.template Solve<MatrixPolicy>(b, x, lower_matrix, upper_matrix);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<double, MatrixPolicy>(x);

  return x.AsVector();
}

// bit to bit variation between CPU and GPU result with randomMatrixVectorOrdering
void verify_gpu_against_cpu()
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
      Group10000CudaDenseMatrix,
      Group10000CudaSparseMatrix,
      micm::CudaLinearSolver<double, Group10000CudaSparseMatrix>>(
      [](const Group10000CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group10000CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group10000CudaSparseMatrix>{ matrix, initial_value };
      },
      10000);

  for (int i = 0; i < cpu_x.size(); i++)
  {
    EXPECT_EQ(cpu_x, gpu_x);
  }
}

TEST(CudaLinearSolver, DenseMatrixVectorOrderingPolicy)
{
  testDenseMatrix<
      Group1CudaDenseMatrix,
      Group1CudaSparseMatrix,
      micm::CudaLinearSolver<double, Group1CudaSparseMatrix, micm::CudaLuDecomposition>>(
      [](const Group1CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1CudaSparseMatrix, micm::CudaLuDecomposition> {
        return micm::CudaLinearSolver<double, Group1CudaSparseMatrix, micm::CudaLuDecomposition>{ matrix, initial_value };
      });
}

TEST(CudaLinearSolver, RandomMatrixVectorOrderingPolicy)
{
  testRandomMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolver<double, Group1CudaSparseMatrix>>(
      [](const Group1CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group1CudaSparseMatrix>{ matrix, initial_value };
      },
      1);
  testRandomMatrix<Group20CudaDenseMatrix, Group20CudaSparseMatrix, micm::CudaLinearSolver<double, Group20CudaSparseMatrix>>(
      [](const Group20CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group20CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group20CudaSparseMatrix>{ matrix, initial_value };
      },
      20);
  testRandomMatrix<Group300CudaDenseMatrix, Group300CudaSparseMatrix, micm::CudaLinearSolver<double, Group300CudaSparseMatrix>>(
      [](const Group300CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group300CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group300CudaSparseMatrix>{ matrix, initial_value };
      },
      300);
  testRandomMatrix<Group4000CudaDenseMatrix, Group4000CudaSparseMatrix, micm::CudaLinearSolver<double, Group4000CudaSparseMatrix>>(
      [](const Group4000CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group4000CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group4000CudaSparseMatrix>{ matrix, initial_value };
      },
      4000);
}

TEST(CudaLinearSolver, DiagonalMatrixVectorOrderingPolicy)
{
  testDiagonalMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolver<double, Group1CudaSparseMatrix>>(
      [](const Group1CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group1CudaSparseMatrix>{ matrix, initial_value };
      },
      1);
  testDiagonalMatrix<Group20CudaDenseMatrix, Group20CudaSparseMatrix, micm::CudaLinearSolver<double, Group20CudaSparseMatrix>>(
      [](const Group20CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group20CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group20CudaSparseMatrix>{ matrix, initial_value };
      },
      20);
  testDiagonalMatrix<Group300CudaDenseMatrix, Group300CudaSparseMatrix, micm::CudaLinearSolver<double, Group300CudaSparseMatrix>>(
      [](const Group300CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group300CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group300CudaSparseMatrix>{ matrix, initial_value };
      },
      300);
  testDiagonalMatrix<Group4000CudaDenseMatrix, Group4000CudaSparseMatrix, micm::CudaLinearSolver<double, Group4000CudaSparseMatrix>>(
      [](const Group4000CudaSparseMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group4000CudaSparseMatrix> {
        return micm::CudaLinearSolver<double, Group4000CudaSparseMatrix>{ matrix, initial_value };
      },
      4000);
}

TEST(CudaLinearSolver, RandomMatrixVectorOrderingForGPU)
{
  verify_gpu_against_cpu();
}
