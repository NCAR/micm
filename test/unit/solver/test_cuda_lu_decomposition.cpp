#pragma once
#include <gtest/gtest.h>

#include <functional>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/cuda_param.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <random>
#include <vector>

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
void gpu_validation(
    const SparseMatrixPolicy<T>& gpu_L,
    const SparseMatrixPolicy<T>& cpu_L,
    const SparseMatrixPolicy<T>& gpu_U,
    const SparseMatrixPolicy<T>& cpu_U)
{
  size_t L_size = cpu_L.AsVector().size();
  size_t U_size = cpu_U.AsVector().size();
  std::vector<T> gpu_L_vector = gpu_L.AsVector();
  std::vector<T> cpu_L_vector = cpu_L.AsVector();
  std::vector<T> gpu_U_vector = gpu_U.AsVector();
  std::vector<T> cpu_U_vector = cpu_U.AsVector();
  for (int i = 0; i < L_size; i++)
  {
    EXPECT_EQ(gpu_L_vector[i], cpu_L_vector[i]);
  };
  for (int j = 0; j < U_size; j++)
  {
    EXPECT_EQ(gpu_U_vector[j], cpu_U_vector[j]);
  };
}

template<template<class> class SparseMatrixPolicy>
void testRandomMatrix(size_t n_grids)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy<double>::create(10).number_of_blocks(n_grids).initial_value(1.0e-30);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.with_element(i, j);

  SparseMatrixPolicy<double> A(builder);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < n_grids; ++i_block)
          A[i_block][i][j] = get_double();

  micm::CudaLuDecomposition gpu_lud(A);
  auto gpu_LU = micm::CudaLuDecomposition::GetLUMatrices(A, 1.0e-30);
  gpu_lud.Decompose<double, SparseMatrixPolicy>(A, gpu_LU.first, gpu_LU.second);
  check_results<double, SparseMatrixPolicy>(
      A, gpu_LU.first, gpu_LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });

  micm::LuDecomposition cpu_lud(A);
  auto cpu_LU = micm::LuDecomposition::GetLUMatrices(A, 1.0e-30);
  cpu_lud.Decompose<double, SparseMatrixPolicy>(A, cpu_LU.first, cpu_LU.second);

  // checking GPU result again CPU
  gpu_validation<double, SparseMatrixPolicy>(gpu_LU.first, cpu_LU.first, gpu_LU.second, cpu_LU.second);
}

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<10>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<100>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1000>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<100000>>;

TEST(CudaLuDecomposition, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix>(10);
  testRandomMatrix<Group2SparseVectorMatrix>(100);
  testRandomMatrix<Group3SparseVectorMatrix>(1000);
  testRandomMatrix<Group4SparseVectorMatrix>(100000);
}