// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

#include <cstdint>
#include <vector>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that calculates the forcing terms on the device
    __global__ void AddForcingTermsKernel(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam forcing_param,
        const ProcessSetParam devstruct)
    {
      // Calculate global thread ID
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      const std::size_t* const __restrict__ d_number_of_reactants = devstruct.number_of_reactants_;
      const std::size_t* __restrict__ d_reactant_ids = devstruct.reactant_ids_;
      const std::size_t* const __restrict__ d_number_of_products = devstruct.number_of_products_;
      const std::size_t* __restrict__ d_product_ids = devstruct.product_ids_;
      const double* __restrict__ d_yields = devstruct.yields_;
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const double* __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* __restrict__ d_state_variables = state_variables_param.d_data_;
      double* __restrict__ d_forcing = forcing_param.d_data_;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t local_tid = tid % cuda_matrix_vector_length;
      const std::size_t group_id = tid / cuda_matrix_vector_length;
      const std::size_t number_of_reactions =
          rate_constants_param.number_of_elements_ / number_of_groups / cuda_matrix_vector_length;

      d_rate_constants += group_id * rate_constants_param.number_of_elements_ / number_of_groups;
      d_state_variables += group_id * state_variables_param.number_of_elements_ / number_of_groups;
      d_forcing += group_id * forcing_param.number_of_elements_ / number_of_groups;

      if (tid < number_of_grid_cells)
      {
        for (int i_rxn = 0; i_rxn < number_of_reactions; ++i_rxn)
        {
          double rate = d_rate_constants[(i_rxn * cuda_matrix_vector_length) + local_tid];
          const std::size_t number_of_reactants = d_number_of_reactants[i_rxn];
          for (int i_react = 0; i_react < number_of_reactants; ++i_react)
          {
            rate *= d_state_variables[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
          }
          for (int i_react = 0; i_react < number_of_reactants; ++i_react)
          {
            d_forcing[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid] -= rate;
          }
          const std::size_t number_of_products = d_number_of_products[i_rxn];
          for (int i_prod = 0; i_prod < number_of_products; ++i_prod)
          {
            d_forcing[d_product_ids[i_prod] * cuda_matrix_vector_length + local_tid] += d_yields[i_prod] * rate;
          }
          d_reactant_ids += d_number_of_reactants[i_rxn];
          d_product_ids += d_number_of_products[i_rxn];
          d_yields += d_number_of_products[i_rxn];
        }  // end of loop over number of reactions
      }  // end of checking a valid CUDA thread id
    }  // end of AddForcingTerms_kernel

    /// This is the CUDA kernel that forms the negative Jacobian matrix (-J) on the device
    __global__ void SubtractJacobianTermsKernel(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      // Calculate global thread ID
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      /// Local device variables
      const ProcessInfoParam* const __restrict__ d_jacobian_process_info = devstruct.jacobian_process_info_;
      const std::size_t* __restrict__ d_reactant_ids = devstruct.jacobian_reactant_ids_;
      const double* __restrict__ d_yields = devstruct.jacobian_yields_;
      const std::size_t* __restrict__ d_jacobian_flat_ids = devstruct.jacobian_flat_ids_;
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;
      const double* __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* __restrict__ d_state_variables = state_variables_param.d_data_;
      double* __restrict__ d_jacobian = jacobian_param.d_data_;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t local_tid = tid % cuda_matrix_vector_length;
      const std::size_t group_id = tid / cuda_matrix_vector_length;

      d_rate_constants += group_id * rate_constants_param.number_of_elements_ / number_of_groups;
      d_state_variables += group_id * state_variables_param.number_of_elements_ / number_of_groups;
      d_jacobian += group_id * jacobian_param.number_of_elements_ / number_of_groups;

      if (tid < number_of_grid_cells)
      {
        // loop over reactions in a grid
        for (int i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
        {
          const ProcessInfoParam& process_info = d_jacobian_process_info[i_proc];
          // Calculate d_rate/d_ind
          double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * cuda_matrix_vector_length) + local_tid];
          for (int i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
          {
            d_rate_d_ind *= d_state_variables[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
          }
          for (int i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
          {
            d_jacobian[*d_jacobian_flat_ids + local_tid] += d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          for (int i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
          {
            std::size_t jacobian_idx = *d_jacobian_flat_ids + local_tid;
            d_jacobian[jacobian_idx] -= d_yields[i_dep] * d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          d_reactant_ids += process_info.number_of_dependent_reactants_;
          d_yields += process_info.number_of_products_;
        }  // end of loop over reactions in a grid cell
      }  // end of checking a CUDA thread id
    }  // end of SubtractJacobianTermsKernel

    /// Optimized CUDA kernel that parallelizes over both grid cells and process_infos (2D grid launch)
    /// blockIdx.x * BLOCK_SIZE + threadIdx.x -> grid cell index (coalesced access)
    /// blockIdx.y -> process_info index (no expensive integer division)
    __global__ void SubtractJacobianTermsKernel2D(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      const std::size_t cell_tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      const std::size_t i_proc = blockIdx.y;

      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      if (cell_tid >= number_of_grid_cells)
        return;

      // Compute group/local indices for vectorized memory layout
      const std::size_t vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + vector_length - 1) / vector_length;
      const std::size_t local_tid = cell_tid % vector_length;
      const std::size_t group_id = cell_tid / vector_length;

      // Base pointers offset by group
      const double* __restrict__ d_rate_constants =
          rate_constants_param.d_data_ + group_id * rate_constants_param.number_of_elements_ / number_of_groups;
      const double* __restrict__ d_state_variables =
          state_variables_param.d_data_ + group_id * state_variables_param.number_of_elements_ / number_of_groups;
      double* __restrict__ d_jacobian =
          jacobian_param.d_data_ + group_id * jacobian_param.number_of_elements_ / number_of_groups;

      // Direct index into process info and offset arrays
      const ProcessInfoParam& process_info = devstruct.jacobian_process_info_[i_proc];
      const std::size_t react_offset = devstruct.jacobian_reactant_offsets_[i_proc];
      const std::size_t flat_offset = devstruct.jacobian_flat_id_offsets_[i_proc];
      const std::size_t yield_offset = devstruct.jacobian_yield_offsets_[i_proc];

      // Compute d_rate/d_ind = rate_constant * product(dependent_reactant_concentrations)
      double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * vector_length) + local_tid];
      for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
      {
        d_rate_d_ind *=
            d_state_variables[(devstruct.jacobian_reactant_ids_[react_offset + i_react] * vector_length) + local_tid];
      }

      // Update jacobian for dependent reactants + independent variable (add d_rate_d_ind)
      for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
      {
        atomicAdd(&d_jacobian[devstruct.jacobian_flat_ids_[flat_offset + i_dep] + local_tid], d_rate_d_ind);
      }

      // Update jacobian for products (subtract yield * d_rate_d_ind)
      const std::size_t product_flat_start = flat_offset + process_info.number_of_dependent_reactants_ + 1;
      for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
      {
        atomicAdd(
            &d_jacobian[devstruct.jacobian_flat_ids_[product_flat_start + i_dep] + local_tid],
            -(devstruct.jacobian_yields_[yield_offset + i_dep] * d_rate_d_ind));
      }
    }  // end of SubtractJacobianTermsKernel2D

    // Kernel-local block size for the compact Jacobian kernel.
    // Smaller than the global BLOCK_SIZE (128) to increase block count and waves/SM,
    // reducing long scoreboard stalls from memory latency.
    static constexpr int JACOBIAN_BLOCK_SIZE = 64;

    /// Optimized kernel using offset-based indexing (no pointer advancement) to break
    /// loop-carried dependencies, enabling instruction-level parallelism across iterations.
    /// Combined with compact types, __ldg(), and JACOBIAN_BLOCK_SIZE=64.
    __global__ void SubtractJacobianTermsKernelCompact(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      std::size_t tid = blockIdx.x * JACOBIAN_BLOCK_SIZE + threadIdx.x;

      const ProcessInfoCompact* const __restrict__ d_process_info = devstruct.jacobian_process_info_compact_;
      const uint16_t* const __restrict__ d_reactant_ids = devstruct.jacobian_reactant_ids_compact_;
      const double* const __restrict__ d_yields = devstruct.jacobian_yields_;
      const uint32_t* const __restrict__ d_flat_ids = devstruct.jacobian_flat_ids_compact_;
      const std::size_t* const __restrict__ d_react_offsets = devstruct.jacobian_reactant_offsets_;
      const std::size_t* const __restrict__ d_flat_offsets = devstruct.jacobian_flat_id_offsets_;
      const std::size_t* const __restrict__ d_yield_offsets = devstruct.jacobian_yield_offsets_;
      const int number_of_grid_cells = static_cast<int>(rate_constants_param.number_of_grid_cells_);
      const int number_of_process_infos = static_cast<int>(devstruct.jacobian_process_info_size_);
      const double* __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* __restrict__ d_state_variables = state_variables_param.d_data_;
      double* __restrict__ d_jacobian = jacobian_param.d_data_;

      const int vector_length = static_cast<int>(state_variables_param.vector_length_);
      const int number_of_groups =
          (number_of_grid_cells + vector_length - 1) / vector_length;
      const int local_tid = static_cast<int>(tid % vector_length);
      const int group_id = static_cast<int>(tid / vector_length);

      d_rate_constants += group_id * rate_constants_param.number_of_elements_ / number_of_groups;
      d_state_variables += group_id * state_variables_param.number_of_elements_ / number_of_groups;
      d_jacobian += group_id * jacobian_param.number_of_elements_ / number_of_groups;

      if (tid < number_of_grid_cells)
      {
        // Offset-based indexing: each iteration independently computes its array positions
        // via pre-computed prefix sums. No loop-carried dependencies between iterations,
        // allowing the compiler to overlap memory loads across iterations (ILP).
        #pragma unroll 4
        for (int i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
        {
          // Read compact process info (8 bytes) via read-only cache
          ProcessInfoCompact pi;
          {
            unsigned long long raw = __ldg(reinterpret_cast<const unsigned long long*>(&d_process_info[i_proc]));
            memcpy(&pi, &raw, sizeof(pi));
          }
          const int num_dep = pi.number_of_dependent_reactants_;
          const int num_prod = pi.number_of_products_;

          // Direct offsets into variable-length arrays (no pointer chasing)
          const int react_off = static_cast<int>(d_react_offsets[i_proc]);
          const int flat_off = static_cast<int>(d_flat_offsets[i_proc]);
          const int yield_off = static_cast<int>(d_yield_offsets[i_proc]);

          // Compute d_rate/d_ind
          double d_rate_d_ind = d_rate_constants[(static_cast<int>(pi.process_id_) * vector_length) + local_tid];
          for (int i_react = 0; i_react < num_dep; ++i_react)
          {
            d_rate_d_ind *= d_state_variables[(__ldg(&d_reactant_ids[react_off + i_react]) * vector_length) + local_tid];
          }

          // Update jacobian for dependent reactants + independent variable
          for (int i_dep = 0; i_dep < num_dep + 1; ++i_dep)
          {
            d_jacobian[__ldg(&d_flat_ids[flat_off + i_dep]) + local_tid] += d_rate_d_ind;
          }

          // Update jacobian for products
          const int prod_flat_off = flat_off + num_dep + 1;
          for (int i_dep = 0; i_dep < num_prod; ++i_dep)
          {
            d_jacobian[__ldg(&d_flat_ids[prod_flat_off + i_dep]) + local_tid] -= __ldg(&d_yields[yield_off + i_dep]) * d_rate_d_ind;
          }
        }
      }
    }  // end of SubtractJacobianTermsKernelCompact

    /// Persistent kernel: each block processes multiple groups sequentially via grid-stride loop.
    /// Block size = vector_length = 128 threads (one thread per cell within a group).
    /// num_blocks is a launch parameter — sweeping it explores the cache-vs-parallelism tradeoff.
    __global__ void SubtractJacobianTermsKernelPersistent(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct,
        int num_groups)
    {
      const ProcessInfoCompact* const __restrict__ d_process_info = devstruct.jacobian_process_info_compact_;
      const uint16_t* const __restrict__ d_reactant_ids = devstruct.jacobian_reactant_ids_compact_;
      const double* const __restrict__ d_yields = devstruct.jacobian_yields_;
      const uint32_t* const __restrict__ d_flat_ids = devstruct.jacobian_flat_ids_compact_;
      const std::size_t* const __restrict__ d_react_offsets = devstruct.jacobian_reactant_offsets_;
      const std::size_t* const __restrict__ d_flat_offsets = devstruct.jacobian_flat_id_offsets_;
      const std::size_t* const __restrict__ d_yield_offsets = devstruct.jacobian_yield_offsets_;
      const int number_of_process_infos = static_cast<int>(devstruct.jacobian_process_info_size_);
      const int vector_length = static_cast<int>(state_variables_param.vector_length_);
      const int local_tid = static_cast<int>(threadIdx.x);

      // Per-group element counts
      const std::size_t rate_const_per_group = rate_constants_param.number_of_elements_ / num_groups;
      const std::size_t state_per_group = state_variables_param.number_of_elements_ / num_groups;
      const std::size_t jac_per_group = jacobian_param.number_of_elements_ / num_groups;

      // Grid-stride loop over groups: each block handles multiple groups sequentially
      for (int group_id = blockIdx.x; group_id < num_groups; group_id += gridDim.x)
      {
        const double* __restrict__ d_rate_constants =
            rate_constants_param.d_data_ + group_id * rate_const_per_group;
        const double* __restrict__ d_state_variables =
            state_variables_param.d_data_ + group_id * state_per_group;
        double* __restrict__ d_jacobian =
            jacobian_param.d_data_ + group_id * jac_per_group;

        // Inner: same as compact kernel body
        for (int i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
        {
          ProcessInfoCompact pi;
          {
            unsigned long long raw = __ldg(reinterpret_cast<const unsigned long long*>(&d_process_info[i_proc]));
            memcpy(&pi, &raw, sizeof(pi));
          }
          const int num_dep = pi.number_of_dependent_reactants_;
          const int num_prod = pi.number_of_products_;

          const int react_off = static_cast<int>(d_react_offsets[i_proc]);
          const int flat_off = static_cast<int>(d_flat_offsets[i_proc]);
          const int yield_off = static_cast<int>(d_yield_offsets[i_proc]);

          double d_rate_d_ind = d_rate_constants[(static_cast<int>(pi.process_id_) * vector_length) + local_tid];
          for (int i_react = 0; i_react < num_dep; ++i_react)
          {
            d_rate_d_ind *= d_state_variables[(__ldg(&d_reactant_ids[react_off + i_react]) * vector_length) + local_tid];
          }

          for (int i_dep = 0; i_dep < num_dep + 1; ++i_dep)
          {
            d_jacobian[__ldg(&d_flat_ids[flat_off + i_dep]) + local_tid] += d_rate_d_ind;
          }

          const int prod_flat_off = flat_off + num_dep + 1;
          for (int i_dep = 0; i_dep < num_prod; ++i_dep)
          {
            d_jacobian[__ldg(&d_flat_ids[prod_flat_off + i_dep]) + local_tid] -= __ldg(&d_yields[yield_off + i_dep]) * d_rate_d_ind;
          }
        }
      }
    }  // end of SubtractJacobianTermsKernelPersistent

    /// This is the function that will copy the constant data
    ///   members of class "CudaProcessSet" to the device,
    ///   except for the "jacobian_flat_id" because it is unknown now
    ProcessSetParam CopyConstData(ProcessSetParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      std::size_t number_of_reactants_bytes = sizeof(std::size_t) * hoststruct.number_of_reactants_size_;
      std::size_t reactant_ids_bytes = sizeof(std::size_t) * hoststruct.reactant_ids_size_;
      std::size_t number_of_products_bytes = sizeof(std::size_t) * hoststruct.number_of_products_size_;
      std::size_t product_ids_bytes = sizeof(std::size_t) * hoststruct.product_ids_size_;
      std::size_t yields_bytes = sizeof(double) * hoststruct.yields_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      ProcessSetParam devstruct;

      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      /// Allocate memory space on the device
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.number_of_reactants_), number_of_reactants_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.reactant_ids_), reactant_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.number_of_products_), number_of_products_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.product_ids_), product_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.yields_), yields_bytes, cuda_stream_id), "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.number_of_reactants_,
              hoststruct.number_of_reactants_,
              number_of_reactants_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.reactant_ids_, hoststruct.reactant_ids_, reactant_ids_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.number_of_products_,
              hoststruct.number_of_products_,
              number_of_products_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.product_ids_, hoststruct.product_ids_, product_ids_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(devstruct.yields_, hoststruct.yields_, yields_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");

      devstruct.number_of_reactants_size_ = hoststruct.number_of_reactants_size_;
      devstruct.reactant_ids_size_ = hoststruct.reactant_ids_size_;
      devstruct.number_of_products_size_ = hoststruct.number_of_products_size_;
      devstruct.product_ids_size_ = hoststruct.product_ids_size_;
      devstruct.yields_size_ = hoststruct.yields_size_;

      return devstruct;
    }

    /// This is the function that will copy the information needed
    /// to calculate Jacobian terms to the device
    void CopyJacobianParams(ProcessSetParam& hoststruct, ProcessSetParam& devstruct)
    {
      /// Calculate the memory space
      std::size_t jacobian_process_info_bytes = sizeof(ProcessInfoParam) * hoststruct.jacobian_process_info_size_;
      std::size_t jacobian_reactant_ids_bytes = sizeof(std::size_t) * hoststruct.jacobian_reactant_ids_size_;
      std::size_t jacobian_product_ids_bytes = sizeof(std::size_t) * hoststruct.jacobian_product_ids_size_;
      std::size_t jacobian_yields_bytes = sizeof(double) * hoststruct.jacobian_yields_size_;
      std::size_t jacobian_flat_ids_bytes = sizeof(std::size_t) * hoststruct.jacobian_flat_ids_size_;

      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      /// Allocate memory space on the device
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_process_info_), jacobian_process_info_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_reactant_ids_), jacobian_reactant_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_product_ids_), jacobian_product_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_yields_), jacobian_yields_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_flat_ids_), jacobian_flat_ids_bytes, cuda_stream_id), "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_process_info_,
              hoststruct.jacobian_process_info_,
              jacobian_process_info_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_reactant_ids_,
              hoststruct.jacobian_reactant_ids_,
              jacobian_reactant_ids_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_product_ids_,
              hoststruct.jacobian_product_ids_,
              jacobian_product_ids_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_yields_,
              hoststruct.jacobian_yields_,
              jacobian_yields_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_flat_ids_,
              hoststruct.jacobian_flat_ids_,
              jacobian_flat_ids_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");

      devstruct.jacobian_process_info_size_ = hoststruct.jacobian_process_info_size_;
      devstruct.jacobian_reactant_ids_size_ = hoststruct.jacobian_reactant_ids_size_;
      devstruct.jacobian_product_ids_size_ = hoststruct.jacobian_product_ids_size_;
      devstruct.jacobian_yields_size_ = hoststruct.jacobian_yields_size_;
      devstruct.jacobian_flat_ids_size_ = hoststruct.jacobian_flat_ids_size_;
    }

    /// Copy pre-computed offset arrays for the optimized 2D kernel to device
    void CopyJacobianOffsets(
        const std::size_t* h_reactant_offsets,
        const std::size_t* h_flat_id_offsets,
        const std::size_t* h_yield_offsets,
        std::size_t offsets_size,
        ProcessSetParam& devstruct)
    {
      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
      std::size_t bytes = sizeof(std::size_t) * offsets_size;

      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_reactant_offsets_), bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_flat_id_offsets_), bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_yield_offsets_), bytes, cuda_stream_id), "cudaMalloc");

      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_reactant_offsets_, h_reactant_offsets, bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_flat_id_offsets_, h_flat_id_offsets, bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_yield_offsets_, h_yield_offsets, bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");

      devstruct.jacobian_offsets_size_ = offsets_size;
    }

    /// Narrow-cast and copy compact Jacobian arrays to device for the bandwidth-optimized kernel
    void CopyJacobianParamsCompact(
        const ProcessInfoParam* h_process_info,
        const std::size_t* h_reactant_ids,
        const std::size_t* h_flat_ids,
        std::size_t process_info_count,
        std::size_t reactant_ids_count,
        std::size_t flat_ids_count,
        ProcessSetParam& devstruct)
    {
      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      // Build compact process info on host
      std::vector<ProcessInfoCompact> compact_info(process_info_count);
      for (std::size_t i = 0; i < process_info_count; ++i)
      {
        compact_info[i].process_id_ = static_cast<uint32_t>(h_process_info[i].process_id_);
        compact_info[i].independent_id_ = static_cast<uint16_t>(h_process_info[i].independent_id_);
        compact_info[i].number_of_dependent_reactants_ = static_cast<uint8_t>(h_process_info[i].number_of_dependent_reactants_);
        compact_info[i].number_of_products_ = static_cast<uint8_t>(h_process_info[i].number_of_products_);
      }
      std::size_t info_bytes = sizeof(ProcessInfoCompact) * process_info_count;
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_process_info_compact_), info_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_process_info_compact_, compact_info.data(), info_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");

      // Narrow-cast flat_ids from size_t to uint32_t
      std::vector<uint32_t> compact_flat_ids(flat_ids_count);
      for (std::size_t i = 0; i < flat_ids_count; ++i)
        compact_flat_ids[i] = static_cast<uint32_t>(h_flat_ids[i]);
      std::size_t flat_bytes = sizeof(uint32_t) * flat_ids_count;
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_flat_ids_compact_), flat_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_flat_ids_compact_, compact_flat_ids.data(), flat_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");

      // Narrow-cast reactant_ids from size_t to uint16_t
      std::vector<uint16_t> compact_reactant_ids(reactant_ids_count);
      for (std::size_t i = 0; i < reactant_ids_count; ++i)
        compact_reactant_ids[i] = static_cast<uint16_t>(h_reactant_ids[i]);
      std::size_t react_bytes = sizeof(uint16_t) * reactant_ids_count;
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_reactant_ids_compact_), react_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_reactant_ids_compact_,
              compact_reactant_ids.data(),
              react_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaProcessSet" on the device
    void FreeConstData(ProcessSetParam& devstruct)
    {
      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      if (devstruct.number_of_reactants_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.number_of_reactants_, cuda_stream_id), "cudaFree");
      if (devstruct.reactant_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.reactant_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.number_of_products_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.number_of_products_, cuda_stream_id), "cudaFree");
      if (devstruct.product_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.product_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.yields_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.yields_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_process_info_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_process_info_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_reactant_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_product_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_product_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_yields_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yields_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_flat_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_reactant_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_flat_id_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_id_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_yield_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yield_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_process_info_compact_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_process_info_compact_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_flat_ids_compact_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_compact_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_reactant_ids_compact_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_compact_, cuda_stream_id), "cudaFree");
    }

    void SubtractJacobianTermsKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SubtractJacobianTermsKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernelDriver

    void AddForcingTermsKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& forcing_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      AddForcingTermsKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, forcing_param, devstruct);
    }  // end of AddForcingTermsKernelDriver

    void SubtractJacobianTermsKernel2DDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      // 2D grid: x-dimension for grid cells, y-dimension for process_infos
      const std::size_t cell_blocks =
          (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      dim3 grid(cell_blocks, devstruct.jacobian_process_info_size_);
      SubtractJacobianTermsKernel2D<<<
          grid,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernel2DDriver

    void SubtractJacobianTermsKernelCompactDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_blocks =
          (rate_constants_param.number_of_grid_cells_ + JACOBIAN_BLOCK_SIZE - 1) / JACOBIAN_BLOCK_SIZE;
      SubtractJacobianTermsKernelCompact<<<
          number_of_blocks,
          JACOBIAN_BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernelCompactDriver

    /// Persistent kernel driver: launches K blocks that grid-stride loop over groups.
    /// Block size = vector_length = 128 (one thread per cell within a group).
    /// num_blocks (K) is parameterized to sweep the cache-vs-parallelism tradeoff.
    void SubtractJacobianTermsKernelPersistentDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct,
        int num_blocks)
    {
      const int vector_length = static_cast<int>(state_variables_param.vector_length_);
      const int num_groups =
          (static_cast<int>(rate_constants_param.number_of_grid_cells_) + vector_length - 1) / vector_length;
      SubtractJacobianTermsKernelPersistent<<<
          num_blocks,
          vector_length,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct, num_groups);
    }  // end of SubtractJacobianTermsKernelPersistentDriver
  }  // namespace cuda
}  // namespace micm
