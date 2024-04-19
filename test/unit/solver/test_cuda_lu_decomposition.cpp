#pragma once
#include <gtest/gtest.h>

#include <functional>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/cuda_param.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/cuda_sparse_matrix.hpp>
#include <random>
#include <vector>

template<typename T, template<class> class SparseMatrixPolicy>
void check_results(
    const SparseMatrixPolicy<T>& A,
    const SparseMatrixPolicy<T>& L,
    const SparseMatrixPolicy<T>& U,
    const std::function<void(const T, const T)> f)
{
  EXPECT_EQ(A.Size(), L.Size());
  EXPECT_EQ(A.Size(), U.Size());
  for (std::size_t i_block = 0; i_block < A.Size(); ++i_block)
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

template<template<class> class CPUSparseMatrixPolicy, template<class> class GPUSparseMatrixPolicy>
void testRandomMatrix(size_t n_grids)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = CPUSparseMatrixPolicy<double>::create(10).number_of_blocks(n_grids).initial_value(1.0e-30);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.with_element(i, j);

  CPUSparseMatrixPolicy<double> cpu_A(builder);
  GPUSparseMatrixPolicy<double> gpu_A(builder);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!cpu_A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < n_grids; ++i_block)
        {
          cpu_A[i_block][i][j] = get_double();
          gpu_A[i_block][i][j] = cpu_A[i_block][i][j];
        }
  
  micm::CudaLuDecomposition gpu_lud(gpu_A);
  auto gpu_LU = micm::CudaLuDecomposition::GetLUMatrices(gpu_A, 1.0e-30);
  gpu_A.CopyToDevice();
  gpu_LU.first.CopyToDevice();
  gpu_LU.second.CopyToDevice();
  gpu_lud.Decompose<double, GPUSparseMatrixPolicy>(gpu_A, gpu_LU.first, gpu_LU.second);
  gpu_LU.first.CopyToHost();
  gpu_LU.second.CopyToHost();
  check_results<double, GPUSparseMatrixPolicy>(
      gpu_A, gpu_LU.first, gpu_LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });

  micm::LuDecomposition cpu_lud = micm::LuDecomposition::Create<double, CPUSparseMatrixPolicy>(cpu_A);
  auto cpu_LU = micm::LuDecomposition::GetLUMatrices<double, CPUSparseMatrixPolicy>(cpu_A, 1.0e-30);
  cpu_lud.Decompose<double, CPUSparseMatrixPolicy>(cpu_A, cpu_LU.first, cpu_LU.second);

  // checking GPU result again CPU
  size_t L_size = cpu_LU.first.AsVector().size();
  size_t U_size = cpu_LU.second.AsVector().size();
  std::vector<double> gpu_L_vector = gpu_LU.first.AsVector();
  std::vector<double> gpu_U_vector = gpu_LU.second.AsVector();
  std::vector<double> cpu_L_vector = cpu_LU.first.AsVector();
  std::vector<double> cpu_U_vector = cpu_LU.second.AsVector();
  for (int i = 0; i < L_size; ++i)
  {
    EXPECT_EQ(gpu_L_vector[i], cpu_L_vector[i]);
  };
  for (int j = 0; j < U_size; ++j)
  {
    EXPECT_EQ(gpu_U_vector[j], cpu_U_vector[j]);
  };
}

template<class T>
using Group1CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group100CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<100>>;
template<class T>
using Group1000CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1000>>;
template<class T>
using Group100000CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<100000>>;

template<class T>
using Group1CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group100CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<100>>;
template<class T>
using Group1000CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<1000>>;
template<class T>
using Group100000CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<100000>>;

TEST(CudaLuDecomposition, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1CPUSparseVectorMatrix, Group1CudaSparseMatrix>(1);
  testRandomMatrix<Group100CPUSparseVectorMatrix, Group100CudaSparseMatrix>(100);
  testRandomMatrix<Group1000CPUSparseVectorMatrix, Group1000CudaSparseMatrix>(1000);
  testRandomMatrix<Group100000CPUSparseVectorMatrix, Group100000CudaSparseMatrix>(100000);
}