#include "../../solver/test_lu_decomposition_policy.hpp"

#include <micm/cuda/solver/cuda_lu_decomposition.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <random>
#include <vector>

template<class CPUSparseMatrixPolicy, class GPUSparseMatrixPolicy>
void testCudaRandomMatrix(size_t n_grids)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = CPUSparseMatrixPolicy::Create(10).SetNumberOfBlocks(n_grids).InitialValue(1.0e-30);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.WithElement(i, j);

  CPUSparseMatrixPolicy cpu_A(builder);
  GPUSparseMatrixPolicy gpu_A(builder);

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
  gpu_lud.Decompose<GPUSparseMatrixPolicy>(gpu_A, gpu_LU.first, gpu_LU.second);
  gpu_LU.first.CopyToHost();
  gpu_LU.second.CopyToHost();
  check_results<typename GPUSparseMatrixPolicy::value_type, GPUSparseMatrixPolicy>(
      gpu_A, gpu_LU.first, gpu_LU.second, [&](const double a, const double b) -> void { EXPECT_LT(std::abs((a-b)/b), 1.0e-07); });

  micm::LuDecomposition cpu_lud = micm::LuDecomposition::Create<CPUSparseMatrixPolicy>(cpu_A);
  auto cpu_LU = micm::LuDecomposition::GetLUMatrices<CPUSparseMatrixPolicy>(cpu_A, 1.0e-30);
  bool singular{ false };
  cpu_lud.Decompose<CPUSparseMatrixPolicy>(cpu_A, cpu_LU.first, cpu_LU.second, singular);

  // checking GPU result again CPU
  size_t L_size = cpu_LU.first.AsVector().size();
  size_t U_size = cpu_LU.second.AsVector().size();
  std::vector<double> gpu_L_vector = gpu_LU.first.AsVector();
  std::vector<double> gpu_U_vector = gpu_LU.second.AsVector();
  std::vector<double> cpu_L_vector = cpu_LU.first.AsVector();
  std::vector<double> cpu_U_vector = cpu_LU.second.AsVector();
  for (int i = 0; i < L_size; ++i)
  {
    auto gpu_L = gpu_L_vector[i];
    auto cpu_L = cpu_L_vector[i];
    EXPECT_LT(std::abs((gpu_L - cpu_L) / cpu_L), 1.0e-10);
  };
  for (int j = 0; j < U_size; ++j)
  {
    auto gpu_U = gpu_U_vector[j];
    auto cpu_U = cpu_U_vector[j];
    EXPECT_LT(std::abs((gpu_U - cpu_U) / cpu_U), 1.0e-10);
  };
}

using Group1CPUSparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group100CPUSparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<100>>;
using Group1000CPUSparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1000>>;
using Group100000CPUSparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<100000>>;

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group100CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<100>>;
using Group1000CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1000>>;
using Group100000CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<100000>>;

TEST(CudaLuDecomposition, RandomMatrixVectorOrdering)
{
  testCudaRandomMatrix<Group1CPUSparseVectorMatrix, Group1CudaSparseMatrix>(1);
  // testCudaRandomMatrix<Group100CPUSparseVectorMatrix, Group100CudaSparseMatrix>(100);
  // testCudaRandomMatrix<Group1000CPUSparseVectorMatrix, Group1000CudaSparseMatrix>(1000);
  // testCudaRandomMatrix<Group100000CPUSparseVectorMatrix, Group100000CudaSparseMatrix>(100000);
}