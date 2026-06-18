// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// CUDA benchmark test for SubtractJacobianTerms using TS1 mechanism const arrays.

#include "subtract_jacobian_terms_ts1_arrays.hpp"

#include <micm/cuda/process/cuda_process_set.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

#include <gtest/gtest.h>

#include <algorithm>
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
static constexpr int MEASURE_ITERATIONS = 1;

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

    // Build prefix-sum offsets needed by SubtractJacobianTermsKernelShared so any
    // thread can jump directly into per-i_proc slices of the const arrays above.
    std::vector<std::size_t> reactant_offsets(JACOBIAN_PROCESS_INFO_SIZE);
    std::vector<std::size_t> yields_offsets(JACOBIAN_PROCESS_INFO_SIZE);
    std::vector<std::size_t> flat_ids_offsets(JACOBIAN_PROCESS_INFO_SIZE);
    std::size_t r_off = 0, y_off = 0, f_off = 0;
    for (std::size_t i = 0; i < JACOBIAN_PROCESS_INFO_SIZE; ++i)
    {
      reactant_offsets[i] = r_off;
      yields_offsets[i] = y_off;
      flat_ids_offsets[i] = f_off;
      r_off += jacobian_process_info_num_dep_reactants[i];
      y_off += jacobian_process_info_num_products[i];
      f_off += jacobian_process_info_num_dep_reactants[i] + 1 + jacobian_process_info_num_products[i];
    }
    devstruct.jacobian_reactant_ids_offsets_ = AllocAndCopy(reactant_offsets.data(), JACOBIAN_PROCESS_INFO_SIZE, stream);
    devstruct.jacobian_yields_offsets_ = AllocAndCopy(yields_offsets.data(), JACOBIAN_PROCESS_INFO_SIZE, stream);
    devstruct.jacobian_flat_ids_offsets_ = AllocAndCopy(flat_ids_offsets.data(), JACOBIAN_PROCESS_INFO_SIZE, stream);

    // Build gather CSR: inverse of jacobian_flat_ids_, keyed by unique Jacobian flat IDs.
    // Each unique flat ID maps to the list of (process_info, sign*yield) contributors.
    // This mirrors ProcessSet::BuildGatherMap() so the same device arrays are available
    // without going through CudaProcessSet.
    struct GatherEntry
    {
      std::size_t jac_flat_id;
      std::size_t proc_idx;
      double coeff;
      std::size_t react_offset;
    };
    std::vector<GatherEntry> gather_entries;
    gather_entries.reserve(JACOBIAN_FLAT_IDS_SIZE);
    {
      std::size_t fp = 0, ro = 0, yo = 0;
      for (std::size_t i = 0; i < JACOBIAN_PROCESS_INFO_SIZE; ++i)
      {
        const std::size_t n_dep = jacobian_process_info_num_dep_reactants[i];
        const std::size_t n_prod = jacobian_process_info_num_products[i];
        for (std::size_t k = 0; k < n_dep + 1; ++k)
          gather_entries.push_back({ jacobian_flat_ids[fp++], i, 1.0, ro });
        for (std::size_t k = 0; k < n_prod; ++k)
          gather_entries.push_back({ jacobian_flat_ids[fp++], i, -jacobian_yields[yo + k], ro });
        ro += n_dep;
        yo += n_prod;
      }
    }
    std::stable_sort(
        gather_entries.begin(),
        gather_entries.end(),
        [](const GatherEntry& a, const GatherEntry& b) { return a.jac_flat_id < b.jac_flat_id; });

    std::vector<std::size_t> h_unique_flat_ids;
    std::vector<std::size_t> h_gather_offsets;
    std::vector<std::size_t> h_proc_idx;
    std::vector<double> h_coeffs;
    std::vector<std::size_t> h_reactant_offset;
    h_gather_offsets.push_back(0);
    h_unique_flat_ids.push_back(gather_entries[0].jac_flat_id);
    for (const auto& e : gather_entries)
    {
      if (e.jac_flat_id != h_unique_flat_ids.back())
      {
        h_unique_flat_ids.push_back(e.jac_flat_id);
        h_gather_offsets.push_back(h_proc_idx.size());
      }
      h_proc_idx.push_back(e.proc_idx);
      h_coeffs.push_back(e.coeff);
      h_reactant_offset.push_back(e.react_offset);
    }
    h_gather_offsets.push_back(h_proc_idx.size());

    devstruct.jac_gather_unique_flat_ids_ = AllocAndCopy(h_unique_flat_ids.data(), h_unique_flat_ids.size(), stream);
    devstruct.jac_gather_offsets_ = AllocAndCopy(h_gather_offsets.data(), h_gather_offsets.size(), stream);
    devstruct.jac_gather_proc_idx_ = AllocAndCopy(h_proc_idx.data(), h_proc_idx.size(), stream);
    devstruct.jac_gather_coeffs_ = AllocAndCopy(h_coeffs.data(), h_coeffs.size(), stream);
    devstruct.jac_gather_reactant_offset_ = AllocAndCopy(h_reactant_offset.data(), h_reactant_offset.size(), stream);
    devstruct.jac_gather_entries_size_ = h_unique_flat_ids.size();

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
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_offsets_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yields_offsets_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_offsets_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_unique_flat_ids_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_offsets_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_proc_idx_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_coeffs_, stream), "cudaFree");
    CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_reactant_offset_, stream), "cudaFree");
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

  // Re-zero after warmup so the measured (and snapshotted) result is a SINGLE solve. The scatter,
  // unrolled, 2D-atomic, and single-pass-gather kernels accumulate with += and would otherwise
  // double-count warmup + measure into the same buffer; the two-pass kernel writes with = (one
  // solve). Zeroing here puts every variant on an identical single-solve basis. It sits outside
  // the timed region, so it does not affect the perf numbers. Assumes MEASURE_ITERATIONS == 1.
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset before measure");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize before measure");

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
  CHECK_CUDA_ERROR(cudaMemcpy(h_jacobian_orig.data(), d_jacobian, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");

  // ---- Benchmark unrolled kernel (#pragma unroll on all loops) ----
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  for (int i = 0; i < WARMUP_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernelUnrolledDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
  }
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after unrolled warmup");

  // Re-zero after warmup so the measured/snapshotted result is a single solve (see note above).
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset before measure");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize before measure");

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

  // ---- Benchmark 2D kernel (grid cells x process_infos, global atomicAdd) ----
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  for (int i = 0; i < WARMUP_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernel2DDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
  }
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after 2D warmup");

  // Re-zero after warmup so the measured/snapshotted result is a single solve (see note above).
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset before measure");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize before measure");

  auto start_2d = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < MEASURE_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsKernel2DDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize during 2D measurement");
  }
  auto end_2d = std::chrono::high_resolution_clock::now();
  double total_2d_ms = std::chrono::duration<double, std::milli>(end_2d - start_2d).count();
  double avg_2d_ms = total_2d_ms / MEASURE_ITERATIONS;

  // Correctness: the 2D kernel uses global atomicAdd to merge contributions from different
  // i_proc threads, so the summation order per Jacobian entry is unspecified. The resulting
  // FP rounding can differ by a few ULPs from the original kernel's serial-loop order.
  // rel_tol = 1e-9 provides ample headroom (same budget used for the former shared kernel).
  std::vector<double> h_jacobian_2d(data.jacobian_total);
  CHECK_CUDA_ERROR(cudaMemcpy(h_jacobian_2d.data(), d_jacobian, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
  std::size_t mismatches_2d = CompareJacobians(h_jacobian_orig, h_jacobian_2d, "2D vs Original", 1e-9);

  // ---- Benchmark single-pass gather kernel (cells x unique Jac entries, no atomics) ----
  // The gather summation order differs from the scatter baseline, so allow 1e-9 tolerance
  // (same as the 2D atomicAdd kernel).
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  for (int i = 0; i < WARMUP_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsGatherKernelDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
  }
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after gather warmup");

  // Re-zero after warmup so the measured/snapshotted result is a single solve (see note above).
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset before measure");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize before measure");

  auto start_gather = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < MEASURE_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsGatherKernelDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize during gather measurement");
  }
  auto end_gather = std::chrono::high_resolution_clock::now();
  double total_gather_ms = std::chrono::duration<double, std::milli>(end_gather - start_gather).count();
  double avg_gather_ms = total_gather_ms / MEASURE_ITERATIONS;

  std::vector<double> h_jacobian_gather(data.jacobian_total);
  CHECK_CUDA_ERROR(
      cudaMemcpy(h_jacobian_gather.data(), d_jacobian, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
  std::size_t mismatches_gather = CompareJacobians(h_jacobian_orig, h_jacobian_gather, "Gather vs Original", 1e-9);

  // ---- Benchmark two-pass gather (Pass1: store d_rate per process; Pass2: gather+accumulate) ----
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  for (int i = 0; i < WARMUP_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsTwoPassGatherDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
  }
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize after two-pass warmup");

  // Re-zero after warmup so the measured/snapshotted result is a single solve (see note above).
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jacobian, 0, jacobian_bytes, stream), "cudaMemset before measure");
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize before measure");

  auto start_twopass = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < MEASURE_ITERATIONS; ++i)
  {
    micm::cuda::SubtractJacobianTermsTwoPassGatherDriver(
        data.rate_constants_param, data.state_variables_param, jacobian_param, data.devstruct);
    CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize during two-pass measurement");
  }
  auto end_twopass = std::chrono::high_resolution_clock::now();
  double total_twopass_ms = std::chrono::duration<double, std::milli>(end_twopass - start_twopass).count();
  double avg_twopass_ms = total_twopass_ms / MEASURE_ITERATIONS;

  std::vector<double> h_jacobian_twopass(data.jacobian_total);
  CHECK_CUDA_ERROR(
      cudaMemcpy(h_jacobian_twopass.data(), d_jacobian, jacobian_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
  std::size_t mismatches_twopass = CompareJacobians(h_jacobian_orig, h_jacobian_twopass, "TwoPass vs Original", 1e-9);

  // ---- Report ----
  const std::size_t n_unique = data.devstruct.jac_gather_entries_size_;
  const std::size_t scratch_mb =
      (NUM_GROUPS * ts1_arrays::JACOBIAN_PROCESS_INFO_SIZE * VECTOR_LENGTH * sizeof(double)) / (1024 * 1024);

  std::cout << "=== SubtractJacobianTerms Performance (TS1) ===" << std::endl;
  std::cout << "  Grid cells:              " << NUM_GRID_CELLS << std::endl;
  std::cout << "  Vector length:           " << VECTOR_LENGTH << std::endl;
  std::cout << "  Num groups:              " << NUM_GROUPS << std::endl;
  std::cout << "  Num reactions:           " << ts1_arrays::NUM_REACTIONS << std::endl;
  std::cout << "  Num species:             " << ts1_arrays::NUM_SPECIES << std::endl;
  std::cout << "  Process info size:       " << ts1_arrays::JACOBIAN_PROCESS_INFO_SIZE << std::endl;
  std::cout << "  Flat IDs size:           " << ts1_arrays::JACOBIAN_FLAT_IDS_SIZE << std::endl;
  std::cout << "  Unique Jac entries:      " << n_unique << std::endl;
  std::cout << "  Two-pass scratch (MB):   " << scratch_mb << std::endl;
  std::cout << "  Warmup iterations:       " << WARMUP_ITERATIONS << std::endl;
  std::cout << "  Measure iterations:      " << MEASURE_ITERATIONS << std::endl;
  std::cout << "  ---" << std::endl;
  std::cout << "  Original (scatter, 1D):  " << avg_orig_ms << " ms/iter" << std::endl;
  std::cout << "  Unrolled (scatter, 1D):  " << avg_unrolled_ms << " ms/iter"
            << "  speedup: " << avg_orig_ms / avg_unrolled_ms << "x" << std::endl;
  std::cout << "  2D (scatter, atomicAdd): " << avg_2d_ms << " ms/iter"
            << "  speedup: " << avg_orig_ms / avg_2d_ms << "x" << std::endl;
  std::cout << "  Gather (single-pass):    " << avg_gather_ms << " ms/iter"
            << "  speedup: " << avg_orig_ms / avg_gather_ms << "x" << std::endl;
  std::cout << "  Gather (two-pass):       " << avg_twopass_ms << " ms/iter"
            << "  speedup: " << avg_orig_ms / avg_twopass_ms << "x" << std::endl;
  std::cout << "===============================================" << std::endl;

  // Cleanup
  CHECK_CUDA_ERROR(cudaFreeAsync(d_jacobian, stream), "cudaFree");
  data.Cleanup(stream);

  EXPECT_EQ(mismatches_unrolled, 0u);
  EXPECT_EQ(mismatches_2d, 0u);
  EXPECT_EQ(mismatches_gather, 0u);
  EXPECT_EQ(mismatches_twopass, 0u);
}

// Dedicated correctness test for both gather kernels using a small problem size so that
// any indexing mistake surfaces without running the full NUM_GRID_CELLS workload.
TEST(CudaSubtractJacobianTerms, GatherKernelsCorrectness)
{
  auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
  using namespace ts1_arrays;

  // Use a single group (128 cells) so discrepancies are easy to trace
  static constexpr std::size_t SMALL_CELLS = VECTOR_LENGTH;
  static constexpr std::size_t SMALL_GROUPS = 1;

  TestDeviceData data;
  data.Setup(stream);

  // Re-create rate_constants / state_variables at the small size (same seed -> same values
  // for the first SMALL_CELLS cells, so we can compare against the large-problem reference
  // by only looking at the first group's worth of Jacobian output).
  // Simpler: run the original kernel on the small problem and use that as reference.
  std::mt19937 gen(99);
  const std::size_t rc_total = SMALL_GROUPS * NUM_REACTIONS * VECTOR_LENGTH;
  const std::size_t sv_total = SMALL_GROUPS * NUM_SPECIES * VECTOR_LENGTH;

  double* d_rc_small = AllocRandomDoubles(rc_total, stream, gen);
  double* d_sv_small = AllocRandomDoubles(sv_total, stream, gen);

  CudaMatrixParam rc_param{};
  rc_param.d_data_ = d_rc_small;
  rc_param.number_of_elements_ = rc_total;
  rc_param.number_of_grid_cells_ = SMALL_CELLS;
  rc_param.vector_length_ = VECTOR_LENGTH;

  CudaMatrixParam sv_param{};
  sv_param.d_data_ = d_sv_small;
  sv_param.number_of_elements_ = sv_total;
  sv_param.number_of_grid_cells_ = SMALL_CELLS;
  sv_param.vector_length_ = VECTOR_LENGTH;

  const std::size_t jac_small_total = SMALL_GROUPS * JACOBIAN_FLAT_BLOCK_SIZE * VECTOR_LENGTH;
  const std::size_t jac_small_bytes = sizeof(double) * jac_small_total;

  double* d_jac_ref = nullptr;
  double* d_jac_gather = nullptr;
  double* d_jac_twopass = nullptr;
  CHECK_CUDA_ERROR(cudaMallocAsync(&d_jac_ref, jac_small_bytes, stream), "cudaMalloc");
  CHECK_CUDA_ERROR(cudaMallocAsync(&d_jac_gather, jac_small_bytes, stream), "cudaMalloc");
  CHECK_CUDA_ERROR(cudaMallocAsync(&d_jac_twopass, jac_small_bytes, stream), "cudaMalloc");

  CudaMatrixParam jac_param{};
  jac_param.number_of_elements_ = jac_small_total;
  jac_param.number_of_grid_cells_ = SMALL_CELLS;
  jac_param.vector_length_ = VECTOR_LENGTH;

  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  // Reference: original scatter kernel
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jac_ref, 0, jac_small_bytes, stream), "cudaMemset");
  jac_param.d_data_ = d_jac_ref;
  micm::cuda::SubtractJacobianTermsKernelDriver(rc_param, sv_param, jac_param, data.devstruct);
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize ref");

  // Single-pass gather
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jac_gather, 0, jac_small_bytes, stream), "cudaMemset");
  jac_param.d_data_ = d_jac_gather;
  micm::cuda::SubtractJacobianTermsGatherKernelDriver(rc_param, sv_param, jac_param, data.devstruct);
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize gather");

  // Two-pass gather
  CHECK_CUDA_ERROR(cudaMemsetAsync(d_jac_twopass, 0, jac_small_bytes, stream), "cudaMemset");
  jac_param.d_data_ = d_jac_twopass;
  micm::cuda::SubtractJacobianTermsTwoPassGatherDriver(rc_param, sv_param, jac_param, data.devstruct);
  CHECK_CUDA_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize twopass");

  // Copy to host and compare
  std::vector<double> h_ref(jac_small_total), h_gather(jac_small_total), h_twopass(jac_small_total);
  CHECK_CUDA_ERROR(cudaMemcpy(h_ref.data(), d_jac_ref, jac_small_bytes, cudaMemcpyDeviceToHost), "D2H ref");
  CHECK_CUDA_ERROR(cudaMemcpy(h_gather.data(), d_jac_gather, jac_small_bytes, cudaMemcpyDeviceToHost), "D2H gather");
  CHECK_CUDA_ERROR(cudaMemcpy(h_twopass.data(), d_jac_twopass, jac_small_bytes, cudaMemcpyDeviceToHost), "D2H twopass");

  std::size_t mm_gather = CompareJacobians(h_ref, h_gather, "Gather (small)", 1e-9);
  std::size_t mm_twopass = CompareJacobians(h_ref, h_twopass, "TwoPass (small)", 1e-9);

  std::cout << "=== GatherKernelsCorrectness (1 group, " << SMALL_CELLS << " cells) ===" << std::endl;
  std::cout << "  Unique Jac entries: " << data.devstruct.jac_gather_entries_size_ << std::endl;
  std::cout << "  Gather mismatches:  " << mm_gather << std::endl;
  std::cout << "  TwoPass mismatches: " << mm_twopass << std::endl;

  CHECK_CUDA_ERROR(cudaFreeAsync(d_rc_small, stream), "cudaFree");
  CHECK_CUDA_ERROR(cudaFreeAsync(d_sv_small, stream), "cudaFree");
  CHECK_CUDA_ERROR(cudaFreeAsync(d_jac_ref, stream), "cudaFree");
  CHECK_CUDA_ERROR(cudaFreeAsync(d_jac_gather, stream), "cudaFree");
  CHECK_CUDA_ERROR(cudaFreeAsync(d_jac_twopass, stream), "cudaFree");
  data.Cleanup(stream);

  EXPECT_EQ(mm_gather, 0u);
  EXPECT_EQ(mm_twopass, 0u);
}
