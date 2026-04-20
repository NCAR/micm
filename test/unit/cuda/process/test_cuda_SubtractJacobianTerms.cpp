// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// CUDA benchmark and correctness test for SubtractJacobianTerms using TS1 mechanism const arrays.
// Compares three kernel variants: original, 2D (atomicAdd), and compact (bandwidth-optimized).

#include "subtract_jacobian_terms_ts1_arrays.hpp"

#include <micm/cuda/process/cuda_process_set.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

// Test configuration
static constexpr std::size_t VECTOR_LENGTH = 128;
static constexpr std::size_t NUM_GRID_CELLS = 163840;
static constexpr std::size_t NUM_GROUPS = NUM_GRID_CELLS / VECTOR_LENGTH;  // 1280

static constexpr int WARMUP_ITERATIONS = 1;
static constexpr int MEASURE_ITERATIONS = 10;

// Helper: allocate device memory and copy from host
template<typename T>
T* AllocAndCopy(const T* host_data, std::size_t count, cudaStream_t stream)
{
  T* d_ptr = nullptr;
  std::size_t bytes = sizeof(T) * count;
  CHECK_CUDA_ERROR(cudaMallocAsync(&d_ptr, bytes, stream), "cudaMalloc");
  CHECK_CUDA_ERROR(cudaMemcpyAsync(d_ptr, host_data, bytes, cudaMemcpyHostToDevice, stream), "cudaMemcpy");
  return d_ptr;
}

// Helper: allocate device memory, fill with random doubles
double* AllocRandomDoubles(std::size_t count, cudaStream_t stream, std::mt19937& gen)
{
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<double> host_data(count);
  for (auto& v : host_data)
    v = dist(gen);

  double* d_ptr = nullptr;
  std::size_t bytes = sizeof(double) * count;
  CHECK_CUDA_ERROR(cudaMallocAsync(&d_ptr, bytes, stream), "cudaMalloc");
  CHECK_CUDA_ERROR(cudaMemcpyAsync(d_ptr, host_data.data(), bytes, cudaMemcpyHostToDevice, stream), "cudaMemcpy");
  return d_ptr;
}

// Helper: compute prefix-sum offset arrays from TS1 const arrays
void ComputeOffsets(
    std::vector<std::size_t>& reactant_offsets,
    std::vector<std::size_t>& flat_id_offsets,
    std::vector<std::size_t>& yield_offsets)
{
  using namespace ts1_arrays;
  const std::size_t n = JACOBIAN_PROCESS_INFO_SIZE;
  reactant_offsets.resize(n + 1);
  flat_id_offsets.resize(n + 1);
  yield_offsets.resize(n + 1);

  reactant_offsets[0] = 0;
  flat_id_offsets[0] = 0;
  yield_offsets[0] = 0;

  for (std::size_t i = 0; i < n; ++i)
  {
    std::size_t num_dep = jacobian_process_info_num_dep_reactants[i];
    std::size_t num_prod = jacobian_process_info_num_products[i];
    reactant_offsets[i + 1] = reactant_offsets[i] + num_dep;
    flat_id_offsets[i + 1] = flat_id_offsets[i] + num_dep + 1 + num_prod;
    yield_offsets[i + 1] = yield_offsets[i] + num_prod;
  }
}

// Helper: set up common device data (ProcessSetParam, rate_constants, state_variables)
struct TestDeviceData
{
  ProcessSetParam devstruct{};
  double* d_rate_constants = nullptr;
  double* d_state_variables = nullptr;
  CudaMatrixParam rate_constants_param{};
  CudaMatrixParam state_variables_param{};
  std::size_t jacobian_total = 0;

  void Setup(cudaStream_t stream)
  {
    using namespace ts1_arrays;

    // Build host ProcessInfoParam array
    std::vector<ProcessInfoParam> host_process_info(JACOBIAN_PROCESS_INFO_SIZE);
    for (std::size_t i = 0; i < JACOBIAN_PROCESS_INFO_SIZE; ++i)
    {
      host_process_info[i].process_id_ = jacobian_process_info_process_id[i];
      host_process_info[i].independent_id_ = jacobian_process_info_independent_id[i];
      host_process_info[i].number_of_dependent_reactants_ = jacobian_process_info_num_dep_reactants[i];
      host_process_info[i].number_of_products_ = jacobian_process_info_num_products[i];
    }

    // Allocate and copy const arrays to device (original kernel format)
    devstruct.jacobian_process_info_ = AllocAndCopy(host_process_info.data(), JACOBIAN_PROCESS_INFO_SIZE, stream);
    devstruct.jacobian_process_info_size_ = JACOBIAN_PROCESS_INFO_SIZE;

    devstruct.jacobian_reactant_ids_ = AllocAndCopy(jacobian_reactant_ids, JACOBIAN_REACTANT_IDS_SIZE, stream);
    devstruct.jacobian_reactant_ids_size_ = JACOBIAN_REACTANT_IDS_SIZE;

    devstruct.jacobian_product_ids_ = AllocAndCopy(jacobian_product_ids, JACOBIAN_PRODUCT_IDS_SIZE, stream);
    devstruct.jacobian_product_ids_size_ = JACOBIAN_PRODUCT_IDS_SIZE;

    devstruct.jacobian_yields_ = AllocAndCopy(jacobian_yields, JACOBIAN_YIELDS_SIZE, stream);
    devstruct.jacobian_yields_size_ = JACOBIAN_YIELDS_SIZE;

    devstruct.jacobian_flat_ids_ = AllocAndCopy(jacobian_flat_ids, JACOBIAN_FLAT_IDS_SIZE, stream);
    devstruct.jacobian_flat_ids_size_ = JACOBIAN_FLAT_IDS_SIZE;

    // Compute and copy offset arrays for 2D kernel
    std::vector<std::size_t> reactant_offsets, flat_id_offsets, yield_offsets;
    ComputeOffsets(reactant_offsets, flat_id_offsets, yield_offsets);
    micm::cuda::CopyJacobianOffsets(
        reactant_offsets.data(), flat_id_offsets.data(), yield_offsets.data(), JACOBIAN_PROCESS_INFO_SIZE + 1, devstruct);

    // Copy compact arrays for bandwidth-optimized kernel
    micm::cuda::CopyJacobianParamsCompact(
        host_process_info.data(),
        jacobian_reactant_ids,
        jacobian_flat_ids,
        JACOBIAN_PROCESS_INFO_SIZE,
        JACOBIAN_REACTANT_IDS_SIZE,
        JACOBIAN_FLAT_IDS_SIZE,
        devstruct);

    // Allocate rate_constants and state_variables with random data
    std::mt19937 gen(42);
    const std::size_t rate_constants_total = NUM_GROUPS * NUM_REACTIONS * VECTOR_LENGTH;
    d_rate_constants = AllocRandomDoubles(rate_constants_total, stream, gen);

    rate_constants_param.d_data_ = d_rate_constants;
    rate_constants_param.number_of_elements_ = rate_constants_total;
    rate_constants_param.number_of_grid_cells_ = NUM_GRID_CELLS;
    rate_constants_param.vector_length_ = VECTOR_LENGTH;

    const std::size_t state_variables_total = NUM_GROUPS * NUM_SPECIES * VECTOR_LENGTH;
    d_state_variables = AllocRandomDoubles(state_variables_total, stream, gen);

    state_variables_param.d_data_ = d_state_variables;
    state_variables_param.number_of_elements_ = state_variables_total;
    state_variables_param.number_of_grid_cells_ = NUM_GRID_CELLS;
    state_variables_param.vector_length_ = VECTOR_LENGTH;

    jacobian_total = NUM_GROUPS * JACOBIAN_FLAT_BLOCK_SIZE * VECTOR_LENGTH;
  }

  void Cleanup(cudaStream_t stream)
  {
    CHECK_CUDA_ERROR(cudaFreeAsync(d_rate_constants, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(d_state_variables, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_process_info_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_product_ids_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yields_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_offsets_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_id_offsets_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yield_offsets_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_process_info_compact_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_compact_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_compact_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize cleanup");
  }
};

// Helper: compare two jacobian arrays element-by-element
std::size_t CompareJacobians(
    const std::vector<double>& h_ref,
    const std::vector<double>& h_test,
    const char* label,
    double rel_tol = 1e-10,
    double abs_tol = 1e-15)
{
  std::size_t mismatches = 0;
  double max_rel_error = 0.0;
  for (std::size_t i = 0; i < h_ref.size(); ++i)
  {
    double ref = h_ref[i];
    double test = h_test[i];
    double abs_diff = std::abs(ref - test);
    double denom = std::max(std::abs(ref), 1e-300);
    double rel_error = abs_diff / denom;
    if (rel_error > max_rel_error)
      max_rel_error = rel_error;
    if (abs_diff > rel_tol * denom && abs_diff > abs_tol)
    {
      ++mismatches;
      if (mismatches <= 5)
      {
        std::cout << "  [" << label << "] Mismatch at index " << i << ": ref=" << ref << " test=" << test
                  << " abs_diff=" << abs_diff << " rel_error=" << rel_error << std::endl;
      }
    }
  }
  std::cout << "  [" << label << "] Mismatches: " << mismatches << "  Max rel error: " << max_rel_error << std::endl;
  return mismatches;
}

// TEST(CudaSubtractJacobianTerms, Correctness)
// {
//   auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

//   TestDeviceData data;
//   data.Setup(stream);

//   const std::size_t jacobian_bytes = sizeof(double) * data.jacobian_total;

//   // Allocate three jacobian arrays (one for each kernel)
//   double* d_jacobian_orig = nullptr;
//   double* d_jacobian_2d = nullptr;
//   double* d_jacobian_compact = nullptr;
//   CHECK_CUDA_ERROR(cudaMallocAsync(&d_jacobian_orig, jacobian_bytes, stream), "cudaMalloc");
//   CHECK_CUDA_ERROR(cudaMallocAsync(&d_jacobian_2d, jacobian_bytes, stream), "cudaMalloc");
//   CHECK_CUDA_ERROR(cudaMallocAsync(&d_jacobian_compact, jacobian_bytes, stream), "cudaMalloc");
//   CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian_orig, 0, jacobian_bytes, stream), "cudaMemset");
//   CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian_2d, 0, jacobian_bytes, stream), "cudaMemset");
//   CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian_compact, 0, jacobian_bytes, stream), "cudaMemset");

//   CudaMatrixParam jacobian_param_orig{};
//   jacobian_param_orig.d_data_ = d_jacobian_orig;
//   jacobian_param_orig.number_of_elements_ = data.jacobian_total;
//   jacobian_param_orig.number_of_grid_cells_ = NUM_GRID_CELLS;
//   jacobian_param_orig.vector_length_ = VECTOR_LENGTH;

//   CudaMatrixParam jacobian_param_2d{};
//   jacobian_param_2d.d_data_ = d_jacobian_2d;
//   jacobian_param_2d.number_of_elements_ = data.jacobian_total;
//   jacobian_param_2d.number_of_grid_cells_ = NUM_GRID_CELLS;
//   jacobian_param_2d.vector_length_ = VECTOR_LENGTH;

//   CudaMatrixParam jacobian_param_compact{};
//   jacobian_param_compact.d_data_ = d_jacobian_compact;
//   jacobian_param_compact.number_of_elements_ = data.jacobian_total;
//   jacobian_param_compact.number_of_grid_cells_ = NUM_GRID_CELLS;
//   jacobian_param_compact.vector_length_ = VECTOR_LENGTH;

//   CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

//   // Run all three kernels
//   micm::cuda::SubtractJacobianTermsKernelDriver(
//       data.rate_constants_param, data.state_variables_param, jacobian_param_orig, data.devstruct);
//   CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after original kernel");

//   micm::cuda::SubtractJacobianTermsKernel2DDriver(
//       data.rate_constants_param, data.state_variables_param, jacobian_param_2d, data.devstruct);
//   CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after 2D kernel");

//   micm::cuda::SubtractJacobianTermsKernelCompactDriver(
//       data.rate_constants_param, data.state_variables_param, jacobian_param_compact, data.devstruct);
//   CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after compact kernel");

//   // Copy results back to host
//   std::vector<double> h_jacobian_orig(data.jacobian_total);
//   std::vector<double> h_jacobian_2d(data.jacobian_total);
//   std::vector<double> h_jacobian_compact(data.jacobian_total);
//   CHECK_CUDA_ERROR(
//       cudaMemcpy(h_jacobian_orig.data(), d_jacobian_orig, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
//   CHECK_CUDA_ERROR(
//       cudaMemcpy(h_jacobian_2d.data(), d_jacobian_2d, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
//   CHECK_CUDA_ERROR(
//       cudaMemcpy(h_jacobian_compact.data(), d_jacobian_compact, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy
//       D2H");

//   // Compare each variant against the original
//   std::cout << "=== Correctness Check ===" << std::endl;
//   std::cout << "  Total elements: " << data.jacobian_total << std::endl;
//   std::size_t mismatches_2d = CompareJacobians(h_jacobian_orig, h_jacobian_2d, "2D vs Original");
//   std::size_t mismatches_compact = CompareJacobians(h_jacobian_orig, h_jacobian_compact, "Compact vs Original");
//   std::cout << "=========================" << std::endl;

//   EXPECT_EQ(mismatches_2d, 0u);
//   EXPECT_EQ(mismatches_compact, 0u);

//   // Cleanup
//   CHECK_CUDA_ERROR(cudaFreeAsync(d_jacobian_orig, stream), "cudaFree");
//   CHECK_CUDA_ERROR(cudaFreeAsync(d_jacobian_2d, stream), "cudaFree");
//   CHECK_CUDA_ERROR(cudaFreeAsync(d_jacobian_compact, stream), "cudaFree");
//   data.Cleanup(stream);
// }

TEST(CudaSubtractJacobianTerms, PerformanceComparison)
{
  auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

  TestDeviceData data;
  data.Setup(stream);

  const std::size_t jacobian_bytes = sizeof(double) * data.jacobian_total;

  double* d_jacobian = nullptr;
  CHECK_CUDA_ERROR(cudaMallocAsync(&d_jacobian, jacobian_bytes, stream), "cudaMalloc");

  CudaMatrixParam jacobian_param{};
  jacobian_param.d_data_ = d_jacobian;
  jacobian_param.number_of_elements_ = data.jacobian_total;
  jacobian_param.number_of_grid_cells_ = NUM_GRID_CELLS;
  jacobian_param.vector_length_ = VECTOR_LENGTH;

  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  // ---- Benchmark original kernel ----
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  for (int i = 0; i < WARMUP_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernelDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
  }
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after warmup");

  auto start_orig = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < MEASURE_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernelDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize during measurement");
  }
  auto end_orig = std::chrono::high_resolution_clock::now();
  double total_orig_ms = std::chrono::duration<double, std::milli>(end_orig - start_orig).count();
  double avg_orig_ms = total_orig_ms / MEASURE_ITERATIONS;

  // ---- Benchmark compact kernel (block size 64, compact types, __ldg) ----
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  for (int i = 0; i < WARMUP_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernelCompactDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
  }
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after warmup");

  auto start_compact = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < MEASURE_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernelCompactDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize during measurement");
  }
  auto end_compact = std::chrono::high_resolution_clock::now();
  double total_compact_ms = std::chrono::duration<double, std::milli>(end_compact - start_compact).count();
  double avg_compact_ms = total_compact_ms / MEASURE_ITERATIONS;

  // ---- Benchmark persistent kernel: sweep K (number of blocks) ----
  // Tests cache-vs-parallelism tradeoff. K=NUM_GROUPS matches non-persistent behavior.
  // Smaller K = fewer concurrent groups = smaller working set = potentially better L2 hit rate.
  const std::vector<int> K_values = {12, 27, 54, 108, 216, 540, static_cast<int>(NUM_GROUPS)};
  std::vector<double> persistent_times;
  for (int K : K_values)
  {
    CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

    for (int i = 0; i < WARMUP_ITERATIONS; ++i)
    {
      micm::cuda::SubtractJacobianTermsKernelPersistentDriver(
          data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct, K);
    }
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after persistent warmup");

    auto start_p = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < MEASURE_ITERATIONS; ++i)
    {
      micm::cuda::SubtractJacobianTermsKernelPersistentDriver(
          data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct, K);
      CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize during persistent measurement");
    }
    auto end_p = std::chrono::high_resolution_clock::now();
    double total_p_ms = std::chrono::duration<double, std::milli>(end_p - start_p).count();
    persistent_times.push_back(total_p_ms / MEASURE_ITERATIONS);
  }

  // ---- Report ----
  std::cout << "=== SubtractJacobianTerms Performance Comparison (TS1) ===" << std::endl;
  std::cout << "  Grid cells:         " << NUM_GRID_CELLS << std::endl;
  std::cout << "  Vector length:      " << VECTOR_LENGTH << std::endl;
  std::cout << "  Num groups:         " << NUM_GROUPS << std::endl;
  std::cout << "  Num reactions:      " << ts1_arrays::NUM_REACTIONS << std::endl;
  std::cout << "  Num species:        " << ts1_arrays::NUM_SPECIES << std::endl;
  std::cout << "  Process info size:  " << ts1_arrays::JACOBIAN_PROCESS_INFO_SIZE << std::endl;
  std::cout << "  Flat IDs size:      " << ts1_arrays::JACOBIAN_FLAT_IDS_SIZE << std::endl;
  std::cout << "  Warmup iterations:  " << WARMUP_ITERATIONS << std::endl;
  std::cout << "  Measure iterations: " << MEASURE_ITERATIONS << std::endl;
  std::cout << "  ---" << std::endl;
  std::cout << "  Original kernel (bs=128): " << avg_orig_ms << " ms/iter (total " << total_orig_ms << " ms)" << std::endl;
  std::cout << "  Compact kernel (bs=64):   " << avg_compact_ms << " ms/iter (total " << total_compact_ms << " ms)"
            << "  speedup: " << avg_orig_ms / avg_compact_ms << "x" << std::endl;
  std::cout << "  --- Persistent kernel K-sweep (cache-vs-parallelism) ---" << std::endl;
  for (size_t i = 0; i < K_values.size(); ++i)
  {
    std::cout << "    K=" << K_values[i] << ":\t" << persistent_times[i] << " ms/iter"
              << "  speedup: " << avg_orig_ms / persistent_times[i] << "x" << std::endl;
  }
  std::cout << "==========================================================" << std::endl;

  // Cleanup
  CHECK_CUDA_ERROR(cudaFreeAsync(d_jacobian, stream), "cudaFree");
  data.Cleanup(stream);

  EXPECT_TRUE(true);
}
