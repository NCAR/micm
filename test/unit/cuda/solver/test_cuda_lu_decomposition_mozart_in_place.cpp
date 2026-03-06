#include "../../solver/test_lu_decomposition_in_place_policy.hpp"

#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.cuh>
#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <random>
#include <vector>

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group100CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<100>>;
using Group1000CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1000>>;
using Group100000CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<100000>>;

constexpr std::size_t number_of_blocks = 163840;
using SparseMatrixPolicy = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<128>>;
using LuDecompositionPolicy = micm::CudaLuDecompositionMozartInPlace;

// TEST(CudaLuDecompositionMozartInPlace, RandomMatrixVectorOrdering)
// {
//   testRandomMatrix<Group1CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1);
//   testRandomMatrix<Group100CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100);
//   testRandomMatrix<Group1000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1000);
//   testRandomMatrix<Group100000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100000);
// }

// TEST(CudaLuDecompositionMozartInPlace, AgnosticToInitialValue)
// {
//   double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
//   for (auto& value : initial_values)
//   {
//     testExtremeValueInitialization<Group1CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1, value);
//     testExtremeValueInitialization<Group100CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100, value);
//     testExtremeValueInitialization<Group1000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1000, value);
//     testExtremeValueInitialization<Group100000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100000, value);
//   }
// }

// TEST(CudaLuDecompositionMozartInPlace, DumpIndexArrays)
// {
//   int n, nnz;

//   FILE* fp =
//       fopen("/glade/derecho/scratch/sunjian/CUDALibrarySamples/cuDSS/simple_uniform_batch/ALU_matrix_output.txt", "r");
//   ASSERT_NE(fp, nullptr) << "Could not open ALU_matrix_output.txt";

//   ASSERT_EQ(fscanf(fp, "%d %d", &nnz, &n), 2);
//   n = n - 1;

//   std::vector<int> csr_offsets(n + 1);
//   std::vector<int> csr_columns(nnz);
//   std::vector<double> csr_values(nnz);

//   for (int i = 0; i < nnz; ++i)
//     ASSERT_EQ(fscanf(fp, "%d", &csr_columns[i]), 1);
//   for (int i = 0; i < n + 1; ++i)
//     ASSERT_EQ(fscanf(fp, "%d", &csr_offsets[i]), 1);
//   for (int i = 0; i < nnz; ++i)
//     ASSERT_EQ(fscanf(fp, "%lf", &csr_values[i]), 1);
//   fclose(fp);

//   auto builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(1).InitialValue(0);
//   for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i)
//     for (std::size_t idx = csr_offsets[i]; idx < static_cast<std::size_t>(csr_offsets[i + 1]); ++idx)
//       builder = builder.WithElement(i, csr_columns[idx]);

//   SparseMatrixPolicy A(builder);

//   LuDecompositionPolicy lud = LuDecompositionPolicy::template Create<SparseMatrixPolicy>(A);

//   // Get the ALU matrix to compute number_of_non_zeros
//   auto ALU = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(A, 0, true);
//   uint32_t number_of_non_zeros = static_cast<uint32_t>(ALU.GroupSize() / SparseMatrixPolicy::GroupVectorSize());

//   const auto& aii_nji_nki = lud.GetAiiNjiNki();
//   const auto& aji = lud.GetAji();
//   const auto& aik_njk = lud.GetAikNjk();
//   const auto& ajk_aji = lud.GetAjkAji();

//   std::size_t nn = aii_nji_nki.size();

//   // Print header with sizes
//   printf(
//       "// n = %zu, aji_size = %zu, aik_njk_size = %zu, ajk_aji_size = %zu, number_of_non_zeros = %u\n",
//       nn,
//       aji.size(),
//       aik_njk.size(),
//       ajk_aji.size(),
//       number_of_non_zeros);

//   // AII
//   printf("static constexpr uint32_t N = %zu;\n", nn);
//   printf("static constexpr uint32_t AJI_SIZE = %zu;\n", aji.size());
//   printf("static constexpr uint32_t AIK_NJK_SIZE = %zu;\n", aik_njk.size());
//   printf("static constexpr uint32_t AJK_AJI_SIZE = %zu;\n", ajk_aji.size());
//   printf("static constexpr uint32_t NUMBER_OF_NON_ZEROS = %u;\n\n", number_of_non_zeros);

//   printf("static constexpr uint32_t AII[%zu] = {", nn);
//   for (std::size_t i = 0; i < nn; ++i)
//     printf("%s%u", i ? ", " : "", static_cast<uint32_t>(std::get<0>(aii_nji_nki[i])));
//   printf("};\n\n");

//   printf("static constexpr uint32_t NJI[%zu] = {", nn);
//   for (std::size_t i = 0; i < nn; ++i)
//     printf("%s%u", i ? ", " : "", static_cast<uint32_t>(std::get<1>(aii_nji_nki[i])));
//   printf("};\n\n");

//   printf("static constexpr uint32_t NKI[%zu] = {", nn);
//   for (std::size_t i = 0; i < nn; ++i)
//     printf("%s%u", i ? ", " : "", static_cast<uint32_t>(std::get<2>(aii_nji_nki[i])));
//   printf("};\n\n");

//   // AJI
//   printf("static constexpr uint32_t AJI[%zu] = {", aji.size());
//   for (std::size_t i = 0; i < aji.size(); ++i)
//     printf("%s%u", i ? ", " : "", static_cast<uint32_t>(aji[i]));
//   printf("};\n\n");

//   // AIK_NJK_PACKED
//   printf("static constexpr uint32_t AIK_NJK_PACKED[%zu] = {", aik_njk.size() * 2);
//   for (std::size_t i = 0; i < aik_njk.size(); ++i)
//   {
//     printf("%s%u, %u", i ? ", " : "", static_cast<uint32_t>(aik_njk[i].first), static_cast<uint32_t>(aik_njk[i].second));
//   }
//   printf("};\n\n");

//   // AJK_AJI_PACKED
//   printf("static constexpr uint32_t AJK_AJI_PACKED[%zu] = {", ajk_aji.size() * 2);
//   for (std::size_t i = 0; i < ajk_aji.size(); ++i)
//   {
//     printf("%s%u, %u", i ? ", " : "", static_cast<uint32_t>(ajk_aji[i].first), static_cast<uint32_t>(ajk_aji[i].second));
//   }
//   printf("};\n\n");
// }

// TEST(CudaLuDecompositionMozartInPlace, TS1_MECHANISM)
// {
//   int n, nnz;

//   // Read matrix data from ALU_matrix_output.txt
//   FILE* fp =
//       fopen("/glade/derecho/scratch/sunjian/CUDALibrarySamples/cuDSS/simple_uniform_batch/ALU_matrix_output.txt", "r");
//   if (!fp)
//   {
//     printf("Error: could not open ALU_matrix_output.txt\n");
//     exit(-1);
//   }

//   if (fscanf(fp, "%d %d", &nnz, &n) != 2)
//   {
//     printf("Error: failed to read matrix sizes from ALU_matrix_output.txt\n");
//     fclose(fp);
//     exit(-1);
//   }
//   n = n - 1;  // actual matrix size is n - 1

//   int* csr_offsets = nullptr;
//   int* csr_columns = nullptr;
//   double* csr_values = nullptr;

//   csr_offsets = (int*)malloc((n + 1) * sizeof(int));
//   csr_columns = (int*)malloc(nnz * sizeof(int));
//   csr_values = (double*)malloc(nnz * sizeof(double));

//   for (int i = 0; i < nnz; ++i)
//   {
//     if (fscanf(fp, "%d", &csr_columns[i]) != 1)
//     {
//       printf("Error: failed to read csr_columns from ALU_matrix_output.txt\n");
//       fclose(fp);
//       free(csr_offsets);
//       free(csr_columns);
//       free(csr_values);
//       exit(-1);
//     }
//   }

//   for (int i = 0; i < n + 1; ++i)
//   {
//     if (fscanf(fp, "%d", &csr_offsets[i]) != 1)
//     {
//       printf("Error: failed to read csr_offsets from ALU_matrix_output.txt\n");
//       fclose(fp);
//       free(csr_offsets);
//       free(csr_columns);
//       free(csr_values);
//       exit(-1);
//     }
//   }

//   for (int i = 0; i < nnz; ++i)
//   {
//     if (fscanf(fp, "%lf", &csr_values[i]) != 1)
//     {
//       printf("Error: failed to read csr_values_h from ALU_matrix_output.txt\n");
//       fclose(fp);
//       free(csr_offsets);
//       free(csr_columns);
//       free(csr_values);
//       exit(-1);
//     }
//   }

//   fclose(fp);

//   auto builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
//   for (std::size_t i = 0; i < n; ++i)
//   {
//     for (std::size_t idx = csr_offsets[i]; idx < csr_offsets[i + 1]; ++idx)
//     {
//       std::size_t j = csr_columns[idx];
//       builder = builder.WithElement(i, j);
//     }
//   }

//   SparseMatrixPolicy A(builder);

//   // Initialize A matrix with values from csr_values
//   for (std::size_t i = 0; i < n; ++i)
//   {
//     for (std::size_t idx = csr_offsets[i]; idx < csr_offsets[i + 1]; ++idx)
//     {
//       std::size_t j = csr_columns[idx];
//       for (std::size_t block = 0; block < number_of_blocks; ++block)
//       {
//         A[block][i][j] = csr_values[j];
//       }
//     }
//   }

//   LuDecompositionPolicy lud = LuDecompositionPolicy::template Create<SparseMatrixPolicy>(A);
//   auto ALU = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(A, 0, false);
//   for (std::size_t i = 0; i < n; ++i)
//     for (std::size_t j = 0; j < n; ++j)
//       if (!A.IsZero(i, j))
//         for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
//           ALU[i_block][i][j] = A[i_block][i][j];

//   CheckCopyToDevice<SparseMatrixPolicy>(ALU);

//   lud.template Decompose<SparseMatrixPolicy>(ALU);

//   CheckCopyToHost<SparseMatrixPolicy>(ALU);
// }

TEST(CudaLuDecompositionMozartInPlace, BenchmarkHardcoded)
{
  int n, nnz;

  FILE* fp =
      fopen("/glade/derecho/scratch/sunjian/CUDALibrarySamples/cuDSS/simple_uniform_batch/ALU_matrix_output.txt", "r");
  ASSERT_NE(fp, nullptr) << "Could not open ALU_matrix_output.txt";

  ASSERT_EQ(fscanf(fp, "%d %d", &nnz, &n), 2);
  n = n - 1;

  std::vector<int> csr_offsets(n + 1);
  std::vector<int> csr_columns(nnz);
  std::vector<double> csr_values(nnz);

  for (int i = 0; i < nnz; ++i)
    ASSERT_EQ(fscanf(fp, "%d", &csr_columns[i]), 1);
  for (int i = 0; i < n + 1; ++i)
    ASSERT_EQ(fscanf(fp, "%d", &csr_offsets[i]), 1);
  for (int i = 0; i < nnz; ++i)
    ASSERT_EQ(fscanf(fp, "%lf", &csr_values[i]), 1);
  fclose(fp);

  auto builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i)
    for (std::size_t idx = csr_offsets[i]; idx < static_cast<std::size_t>(csr_offsets[i + 1]); ++idx)
      builder = builder.WithElement(i, csr_columns[idx]);

  SparseMatrixPolicy A(builder);

  for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i)
    for (std::size_t idx = csr_offsets[i]; idx < static_cast<std::size_t>(csr_offsets[i + 1]); ++idx)
      for (std::size_t block = 0; block < number_of_blocks; ++block)
        A[block][i][csr_columns[idx]] = csr_values[csr_columns[idx]];

  LuDecompositionPolicy lud = LuDecompositionPolicy::template Create<SparseMatrixPolicy>(A);
  auto ALU_orig = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(A, 0, false);
  uint32_t nnz_per_block = static_cast<uint32_t>(ALU_orig.GroupSize() / SparseMatrixPolicy::GroupVectorSize());

  for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i)
    for (std::size_t j = 0; j < static_cast<std::size_t>(n); ++j)
      if (!A.IsZero(i, j))
        for (std::size_t blk = 0; blk < number_of_blocks; ++blk)
          ALU_orig[blk][i][j] = A[blk][i][j];

  CheckCopyToDevice<SparseMatrixPolicy>(ALU_orig);

  // Make a copy for the hardcoded kernel
  auto ALU_hard = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(A, 0, false);
  for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i)
    for (std::size_t j = 0; j < static_cast<std::size_t>(n); ++j)
      if (!A.IsZero(i, j))
        for (std::size_t blk = 0; blk < number_of_blocks; ++blk)
          ALU_hard[blk][i][j] = A[blk][i][j];

  CheckCopyToDevice<SparseMatrixPolicy>(ALU_hard);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  const int warmup = 1;
  const int iterations = 10;

  // Benchmark original kernel
  for (int i = 0; i < warmup; ++i)
    lud.template Decompose<SparseMatrixPolicy>(ALU_orig);
  cudaDeviceSynchronize();

  cudaEventRecord(start);
  for (int i = 0; i < iterations; ++i)
    lud.template Decompose<SparseMatrixPolicy>(ALU_orig);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float ms_original = 0;
  cudaEventElapsedTime(&ms_original, start, stop);

  // Benchmark hardcoded kernel
  auto ALU_hard_param = ALU_hard.AsDeviceParam();
  for (int i = 0; i < warmup; ++i)
    micm::cuda::DecomposeKernelHardcodedDriver(ALU_hard_param, nnz_per_block);
  cudaDeviceSynchronize();

  cudaEventRecord(start);
  for (int i = 0; i < iterations; ++i)
    micm::cuda::DecomposeKernelHardcodedDriver(ALU_hard_param, nnz_per_block);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float ms_hardcoded = 0;
  cudaEventElapsedTime(&ms_hardcoded, start, stop);

  printf("Original kernel:  %.3f ms (%d iterations)\n", ms_original, iterations);
  printf("Hardcoded kernel: %.3f ms (%d iterations)\n", ms_hardcoded, iterations);
  printf("Speedup: %.2fx\n", ms_original / ms_hardcoded);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}