#include <gtest/gtest.h>

#include <functional>
#include <micm/solver/linear_solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/error_policies.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>

template<typename T, template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void check_results(
    const SparseMatrixPolicy<T> A,
    const MatrixPolicy<T> b,
    const MatrixPolicy<T> x,
    const std::function<void(const T, const T)> f)
{
  T result;
  EXPECT_EQ(A.size(), b.size());
  EXPECT_EQ(A.size(), x.size());
  for (std::size_t i_block = 0; i_block < A.size(); ++i_block)
  {
    for (std::size_t i = 0; i < A[i_block].size(); ++i)
    {
      result = 0.0;
      for (std::size_t j = 0; j < A[i_block].size(); ++j)
        if (!A.IsZero(i, j))
          result += A[i_block][i][j] * x[i_block][j];
      f(b[i_block][i], result);
    }
  }
}

template<typename T, template<class> class SparseMatrixPolicy>
void print_matrix(const SparseMatrixPolicy<T>& matrix, std::size_t width)
{
  for (std::size_t i_block = 0; i_block < matrix.size(); ++i_block)
  {
    std::cout << "block: " << i_block << std::endl;
    for (std::size_t i = 0; i < matrix[i_block].size(); ++i)
    {
      for (std::size_t j = 0; j < matrix[i_block][i].size(); ++j)
      {
        if (matrix.IsZero(i, j))
        {
          std::cout << " " << std::setfill('-') << std::setw(width) << "-"
                    << " ";
        }
        else
        {
          std::cout << " " << std::setfill(' ') << std::setw(width) << matrix[i_block][i][j] << " ";
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void testDenseMatrix()
{
  SparseMatrixPolicy<double> A = SparseMatrixPolicy<double>(SparseMatrixPolicy<double>::create(3)
                                                                .initial_value(1.0e-30)
                                                                .with_element(0, 0)
                                                                .with_element(0, 1)
                                                                .with_element(0, 2)
                                                                .with_element(1, 0)
                                                                .with_element(1, 1)
                                                                .with_element(1, 2)
                                                                .with_element(2, 0)
                                                                .with_element(2, 1)
                                                                .with_element(2, 2));
  MatrixPolicy<double> b(1, 3, 0.0);
  MatrixPolicy<double> x(1, 3, 100.0);

  A[0][0][0] = 2;
  A[0][0][1] = -1;
  A[0][0][2] = -2;
  A[0][1][0] = -4;
  A[0][1][1] = 6;
  A[0][1][2] = 3;
  A[0][2][0] = -4;
  A[0][2][1] = -2;
  A[0][2][2] = 8;

  b[0][0] = 23;
  b[0][1] = 42;
  b[0][2] = 9;

  micm::LinearSolver<double, SparseMatrixPolicy> solver(A, 1.0e-30);
  solver.Factor(A);
  solver.template Solve<MatrixPolicy>(b, x);
  check_results<double, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void testRandomMatrix()
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy<double>::create(10).number_of_blocks(5).initial_value(1.0e-30);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.with_element(i, j);

  SparseMatrixPolicy<double> A(builder);
  MatrixPolicy<double> b(5, 10, 0.0);
  MatrixPolicy<double> x(5, 10, 100.0);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < 5; ++i_block)
          A[i_block][i][j] = get_double();

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t i_block = 0; i_block < 5; ++i_block)
      b[i_block][i] = get_double();

  micm::LinearSolver<double, SparseMatrixPolicy> solver(A, 1.0e-30);
  solver.Factor(A);
  solver.template Solve<MatrixPolicy>(b, x);
  check_results<double, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void testDiagonalMatrix()
{
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy<double>::create(6).number_of_blocks(5).initial_value(1.0e-30);
  for (std::size_t i = 0; i < 6; ++i)
    builder = builder.with_element(i, i);

  SparseMatrixPolicy<double> A(builder);
  MatrixPolicy<double> b(5, 6, 0.0);
  MatrixPolicy<double> x(5, 6, 100.0);

  for (std::size_t i = 0; i < 6; ++i)
    for (std::size_t i_block = 0; i_block < 5; ++i_block)
      A[i_block][i][i] = get_double();

  micm::LinearSolver<double, SparseMatrixPolicy> solver(A, 1.0e-30);
  solver.Factor(A);
  solver.template Solve<MatrixPolicy>(b, x);
  check_results<double, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<class T>
using SparseMatrix = micm::SparseMatrix<T>;
TEST(LinearSolver, DenseMatrixStandardOrdering)
{
  testDenseMatrix<micm::Matrix, SparseMatrix>();
}

TEST(LinearSolver, RandomMatrixStandardOrdering)
{
  testRandomMatrix<micm::Matrix, SparseMatrix>();
}

TEST(LinearSolver, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<micm::Matrix, SparseMatrix>();
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
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::InvalidArgumentPolicy, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::InvalidArgumentPolicy, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::InvalidArgumentPolicy, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::InvalidArgumentPolicy, micm::SparseMatrixVectorOrdering<4>>;

TEST(LinearSolver, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1VectorMatrix, Group1SparseVectorMatrix>();
  testDenseMatrix<Group2VectorMatrix, Group2SparseVectorMatrix>();
  testDenseMatrix<Group3VectorMatrix, Group3SparseVectorMatrix>();
  testDenseMatrix<Group4VectorMatrix, Group4SparseVectorMatrix>();
}

TEST(LinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix>();
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix>();
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix>();
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix>();
}

TEST(LinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix>();
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix>();
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix>();
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix>();
}
