#include <gtest/gtest.h>

#include <functional>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <random>

template<typename T, template<class> class SparseMatrixPolicy>
void check_results(
    const SparseMatrixPolicy<T>& A,
    const SparseMatrixPolicy<T>& L,
    const SparseMatrixPolicy<T>& U,
    const std::function<void(const T, const T)> f)
{
  EXPECT_EQ(A.size(), L.size());
  EXPECT_EQ(A.size(), U.size());
  for (std::size_t i_block = 0; i_block < A.size(); ++i_block)
  {
    for (std::size_t i = 0; i < A[i_block].size(); ++i)
    {
      for (std::size_t j = 0; j < A[i_block].size(); ++j)
      {
        T result{};
        for (std::size_t k = 0; k < A[i_block].size(); ++k)
        {
          if (!(L.IsZero(i, k) || U.IsZero(k, j)))
          {
            result += L[i_block][i][k] * U[i_block][k][j];
          }
        }
        // Make sure these are actually triangular matrices
        EXPECT_TRUE(i >= j || L.IsZero(i, j));
        EXPECT_TRUE(j >= i || U.IsZero(i, j));
        if (A.IsZero(i, j))
        {
          f(result, T{});
        }
        else
        {
          f(result, A[i_block][i][j]);
        }
      }
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

// tests example from https://www.geeksforgeeks.org/doolittle-algorithm-lu-decomposition/
template<template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
void testDenseMatrix(const std::function<LuDecompositionPolicy(const SparseMatrixPolicy<double>&)> create_lu_decomp)
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

  A[0][0][0] = 2;
  A[0][0][1] = -1;
  A[0][0][2] = -2;
  A[0][1][0] = -4;
  A[0][1][1] = 6;
  A[0][1][2] = 3;
  A[0][2][0] = -4;
  A[0][2][1] = -2;
  A[0][2][2] = 8;

  LuDecompositionPolicy lud = create_lu_decomp(A);
  auto LU = micm::LuDecomposition::GetLUMatrices(A, 1.0e-30);
  lud.template Decompose<double, SparseMatrixPolicy>(A, LU.first, LU.second);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
void testSingularMatrix(const std::function<LuDecompositionPolicy(const SparseMatrixPolicy<double>&)> create_lu_decomp)
{
  SparseMatrixPolicy<double> A = SparseMatrixPolicy<double>(SparseMatrixPolicy<double>::create(2)
                                                                .initial_value(1.0e-30)
                                                                .with_element(0, 0)
                                                                .with_element(0, 1)
                                                                .with_element(1, 0)
                                                                .with_element(1, 1));

  A[0][0][0] = 0;
  A[0][0][1] = 1;
  A[0][1][0] = 1;
  A[0][1][1] = 1;

  LuDecompositionPolicy lud = create_lu_decomp(A);
  auto LU = micm::LuDecomposition::GetLUMatrices(A, 1.0E-30);
  bool is_singular{ false };
  lud.template Decompose<double, SparseMatrixPolicy>(A, LU.first, LU.second, is_singular);
  EXPECT_TRUE(is_singular);
  A[0][0][0] = 12;
  lud.template Decompose<double, SparseMatrixPolicy>(A, LU.first, LU.second, is_singular);
  EXPECT_FALSE(is_singular);
}

template<template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
void testRandomMatrix(
    const std::function<LuDecompositionPolicy(const SparseMatrixPolicy<double>&)> create_lu_decomp,
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

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = get_double();

  LuDecompositionPolicy lud = create_lu_decomp(A);
  auto LU = micm::LuDecomposition::GetLUMatrices(A, 1.0e-30);
  lud.template Decompose<double, SparseMatrixPolicy>(A, LU.first, LU.second);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
void testDiagonalMatrix(
    const std::function<LuDecompositionPolicy(const SparseMatrixPolicy<double>&)> create_lu_decomp,
    std::size_t number_of_blocks)
{
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy<double>::create(6).number_of_blocks(number_of_blocks).initial_value(1.0e-30);
  for (std::size_t i = 0; i < 6; ++i)
    builder = builder.with_element(i, i);

  SparseMatrixPolicy<double> A(builder);

  for (std::size_t i = 0; i < 6; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      A[i_block][i][i] = get_double();

  LuDecompositionPolicy lud = create_lu_decomp(A);
  auto LU = micm::LuDecomposition::GetLUMatrices(A, 1.0e-30);
  lud.template Decompose<double, SparseMatrixPolicy>(A, LU.first, LU.second);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}