// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// CUDA benchmark test for SubtractJacobianTerms using TS1 mechanism const arrays.

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

  // Snapshot the original kernel's output for correctness comparison below
  std::vector<double> h_jacobian_orig(data.jacobian_total);
  CHECK_CUDA_ERROR(
      cudaMemcpy(h_jacobian_orig.data(), d_jacobian, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");

  // ---- Benchmark unrolled kernel (#pragma unroll on all loops) ----
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  for (int i = 0; i < WARMUP_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernelUnrolledDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
  }
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after unrolled warmup");

  auto start_unrolled = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < MEASURE_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernelUnrolledDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize during unrolled measurement");
  }
  auto end_unrolled = std::chrono::high_resolution_clock::now();
  double total_unrolled_ms = std::chrono::duration<double, std::milli>(end_unrolled - start_unrolled).count();
  double avg_unrolled_ms = total_unrolled_ms / MEASURE_ITERATIONS;

  // Correctness check: unrolled output should match the original element-for-element
  std::vector<double> h_jacobian_unrolled(data.jacobian_total);
  CHECK_CUDA_ERROR(
      cudaMemcpy(h_jacobian_unrolled.data(), d_jacobian, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
  std::size_t mismatches_unrolled = CompareJacobians(h_jacobian_orig, h_jacobian_unrolled, "Unrolled vs Original");

  // ---- Report ----
  std::cout << "=== SubtractJacobianTerms Performance (TS1) ===" << std::endl;
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
  std::cout << "  Unrolled kernel (bs=128): " << avg_unrolled_ms << " ms/iter (total " << total_unrolled_ms << " ms)"
            << "  speedup: " << avg_orig_ms / avg_unrolled_ms << "x" << std::endl;
  std::cout << "===============================================" << std::endl;

  // Cleanup
  CHECK_CUDA_ERROR(cudaFreeAsync(d_jacobian, stream), "cudaFree");
  data.Cleanup(stream);

  EXPECT_EQ(mismatches_unrolled, 0u);
}
