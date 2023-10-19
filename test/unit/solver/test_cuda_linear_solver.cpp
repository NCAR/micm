#pragma once
#include <gtest/gtest.h>
#include <functional>
#include <random>
#include <micm/solver/linear_solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <micm/util/cuda_param.hpp>
#include <micm/solver/cuda_linear_solver.hpp>

#include "test_linear_solver_policy.hpp"
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

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testRandomMatrix(
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
  solver.Factor(A);
  solver.template Solve<MatrixPolicy>(b, x);
  check_results<double, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}


