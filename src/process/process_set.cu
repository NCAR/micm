// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

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

    /// Unrolled variant of SubtractJacobianTermsKernel. Same data layout and access pattern
    /// as the original (pointer advancement over the const arrays); only difference is
    /// #pragma unroll hints on every loop.
    __global__ void SubtractJacobianTermsKernelUnrolled(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

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
        for (int i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
        {
          const ProcessInfoParam& process_info = d_jacobian_process_info[i_proc];
          double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * cuda_matrix_vector_length) + local_tid];
          #pragma unroll
          for (int i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
          {
            d_rate_d_ind *= d_state_variables[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
          }
          #pragma unroll
          for (int i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
          {
            d_jacobian[*d_jacobian_flat_ids + local_tid] += d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          #pragma unroll
          for (int i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
          {
            std::size_t jacobian_idx = *d_jacobian_flat_ids + local_tid;
            d_jacobian[jacobian_idx] -= d_yields[i_dep] * d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          d_reactant_ids += process_info.number_of_dependent_reactants_;
          d_yields += process_info.number_of_products_;
        }
      }
    }  // end of SubtractJacobianTermsKernelUnrolled

    /// 2D kernel that parallelizes both grid cells (X) and process_infos (Y) independently.
    /// Grid: (ceil(num_cells/BLOCK_DIM_X), ceil(num_proc_infos/BLOCK_DIM_Y)).
    /// Block: (BLOCK_DIM_X, BLOCK_DIM_Y) with BLOCK_DIM_X * BLOCK_DIM_Y == BLOCK_SIZE.
    /// Each thread handles exactly one (cell_id, i_proc) pair. Contributions from different
    /// i_proc threads targeting the same Jacobian entry are merged via global atomicAdd.
    /// Requires the three jacobian_*_offsets_ arrays in devstruct to be populated so each
    /// thread can locate its i_proc slice in O(1).
    __global__ void SubtractJacobianTermsKernel2D(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      const std::size_t cell_id = blockIdx.x * blockDim.x + threadIdx.x;
      const std::size_t i_proc = blockIdx.y * blockDim.y + threadIdx.y;

      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;

      if (cell_id >= number_of_grid_cells || i_proc >= number_of_process_infos)
        return;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t group_id = cell_id / cuda_matrix_vector_length;
      const std::size_t local_tid = cell_id % cuda_matrix_vector_length;

      const double* __restrict__ d_rate_constants =
          rate_constants_param.d_data_ + group_id * (rate_constants_param.number_of_elements_ / number_of_groups);
      const double* __restrict__ d_state_variables =
          state_variables_param.d_data_ + group_id * (state_variables_param.number_of_elements_ / number_of_groups);
      double* __restrict__ d_jacobian =
          jacobian_param.d_data_ + group_id * (jacobian_param.number_of_elements_ / number_of_groups);

      const ProcessInfoParam process_info = devstruct.jacobian_process_info_[i_proc];
      const std::size_t* __restrict__ my_reactant_ids =
          devstruct.jacobian_reactant_ids_ + devstruct.jacobian_reactant_ids_offsets_[i_proc];
      const double* __restrict__ my_yields = devstruct.jacobian_yields_ + devstruct.jacobian_yields_offsets_[i_proc];
      const std::size_t* __restrict__ my_flat_ids =
          devstruct.jacobian_flat_ids_ + devstruct.jacobian_flat_ids_offsets_[i_proc];

      double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * cuda_matrix_vector_length) + local_tid];
      for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
      {
        d_rate_d_ind *= d_state_variables[(my_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
      }

      const std::size_t n_pos = process_info.number_of_dependent_reactants_ + 1;
      for (std::size_t i_dep = 0; i_dep < n_pos; ++i_dep)
      {
        atomicAdd(&d_jacobian[my_flat_ids[i_dep] + local_tid], d_rate_d_ind);
      }
      for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
      {
        atomicAdd(&d_jacobian[my_flat_ids[n_pos + i_dep] + local_tid], -my_yields[i_dep] * d_rate_d_ind);
      }
    }  // end of SubtractJacobianTermsKernel2D

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

      /// Gather CSR arrays
      std::size_t jac_gather_unique_flat_ids_bytes = sizeof(std::size_t) * hoststruct.jac_gather_entries_size_;
      std::size_t jac_gather_offsets_bytes = sizeof(std::size_t) * (hoststruct.jac_gather_entries_size_ + 1);
      std::size_t jac_gather_per_source_bytes = sizeof(std::size_t) * hoststruct.jacobian_flat_ids_size_;
      std::size_t jac_gather_coeffs_bytes = sizeof(double) * hoststruct.jacobian_flat_ids_size_;

      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jac_gather_unique_flat_ids_), jac_gather_unique_flat_ids_bytes, cuda_stream_id),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jac_gather_offsets_), jac_gather_offsets_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jac_gather_proc_idx_), jac_gather_per_source_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jac_gather_coeffs_), jac_gather_coeffs_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jac_gather_reactant_offset_), jac_gather_per_source_bytes, cuda_stream_id),
          "cudaMalloc");

      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jac_gather_unique_flat_ids_,
              hoststruct.jac_gather_unique_flat_ids_,
              jac_gather_unique_flat_ids_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jac_gather_offsets_,
              hoststruct.jac_gather_offsets_,
              jac_gather_offsets_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jac_gather_proc_idx_,
              hoststruct.jac_gather_proc_idx_,
              jac_gather_per_source_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jac_gather_coeffs_,
              hoststruct.jac_gather_coeffs_,
              jac_gather_coeffs_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jac_gather_reactant_offset_,
              hoststruct.jac_gather_reactant_offset_,
              jac_gather_per_source_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");

      devstruct.jac_gather_entries_size_ = hoststruct.jac_gather_entries_size_;
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
      if (devstruct.jacobian_reactant_ids_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_yields_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yields_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_flat_ids_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jac_gather_unique_flat_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_unique_flat_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.jac_gather_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jac_gather_proc_idx_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_proc_idx_, cuda_stream_id), "cudaFree");
      if (devstruct.jac_gather_coeffs_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_coeffs_, cuda_stream_id), "cudaFree");
      if (devstruct.jac_gather_reactant_offset_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jac_gather_reactant_offset_, cuda_stream_id), "cudaFree");
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

    void SubtractJacobianTermsKernelUnrolledDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SubtractJacobianTermsKernelUnrolled<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernelUnrolledDriver

    void SubtractJacobianTermsKernel2DDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      constexpr std::size_t BLOCK_DIM_X = 32;
      constexpr std::size_t BLOCK_DIM_Y = 4;
      static_assert(BLOCK_DIM_X * BLOCK_DIM_Y == BLOCK_SIZE, "2D block must equal BLOCK_SIZE");

      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;

      const dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
      const dim3 grid(
          (number_of_grid_cells + BLOCK_DIM_X - 1) / BLOCK_DIM_X,
          (number_of_process_infos + BLOCK_DIM_Y - 1) / BLOCK_DIM_Y);

      SubtractJacobianTermsKernel2D<<<grid, block, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernel2DDriver

    /// Single-pass gather kernel: parallelizes on (grid cell, unique Jacobian entry).
    /// Each thread owns exactly one Jacobian entry and loops over the CSR contributors to
    /// accumulate the total. Single write per entry — no atomics needed.
    /// Trade-off: recomputes d_rate_d_ind for each contributing process per entry (duplicate
    /// compute when multiple entries share the same process).
    __global__ void SubtractJacobianTermsGatherKernel(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      const std::size_t cell_id = blockIdx.x * blockDim.x + threadIdx.x;
      const std::size_t jac_idx = blockIdx.y * blockDim.y + threadIdx.y;

      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      if (cell_id >= number_of_grid_cells || jac_idx >= devstruct.jac_gather_entries_size_)
        return;

      const std::size_t VL = state_variables_param.vector_length_;
      const std::size_t number_of_groups = (number_of_grid_cells + VL - 1) / VL;
      const std::size_t group_id = cell_id / VL;
      const std::size_t local_tid = cell_id % VL;

      const double* __restrict__ d_rate_constants =
          rate_constants_param.d_data_ + group_id * (rate_constants_param.number_of_elements_ / number_of_groups);
      const double* __restrict__ d_state_variables =
          state_variables_param.d_data_ + group_id * (state_variables_param.number_of_elements_ / number_of_groups);
      double* __restrict__ d_jacobian =
          jacobian_param.d_data_ + group_id * (jacobian_param.number_of_elements_ / number_of_groups);

      const std::size_t start = devstruct.jac_gather_offsets_[jac_idx];
      const std::size_t end = devstruct.jac_gather_offsets_[jac_idx + 1];

      double total = 0.0;
      for (std::size_t s = start; s < end; ++s)
      {
        const std::size_t i_proc = devstruct.jac_gather_proc_idx_[s];
        const std::size_t react_off = devstruct.jac_gather_reactant_offset_[s];
        const ProcessInfoParam pi = devstruct.jacobian_process_info_[i_proc];
        double d_rate = d_rate_constants[pi.process_id_ * VL + local_tid];
        for (std::size_t r = 0; r < pi.number_of_dependent_reactants_; ++r)
          d_rate *= d_state_variables[devstruct.jacobian_reactant_ids_[react_off + r] * VL + local_tid];
        total += devstruct.jac_gather_coeffs_[s] * d_rate;
      }
      // Single write per (jac_idx, local_tid) — no atomic needed
      d_jacobian[devstruct.jac_gather_unique_flat_ids_[jac_idx] + local_tid] += total;
    }  // end of SubtractJacobianTermsGatherKernel

    void SubtractJacobianTermsGatherKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      constexpr std::size_t BLOCK_DIM_X = 32;
      constexpr std::size_t BLOCK_DIM_Y = 4;
      static_assert(BLOCK_DIM_X * BLOCK_DIM_Y == BLOCK_SIZE, "2D block must equal BLOCK_SIZE");

      const std::size_t N_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t N_unique = devstruct.jac_gather_entries_size_;

      const dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
      const dim3 grid((N_cells + BLOCK_DIM_X - 1) / BLOCK_DIM_X, (N_unique + BLOCK_DIM_Y - 1) / BLOCK_DIM_Y);

      SubtractJacobianTermsGatherKernel<<<
          grid,
          block,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsGatherKernelDriver

    /// Two-pass gather, Pass 1: compute d_rate_d_ind for every process_info and store in
    /// scratch. One value per (group, process, local_tid). No writes to the Jacobian.
    /// Scratch layout: scratch[group_id * N_proc * VL + i_proc * VL + local_tid].
    __global__ void SubtractJacobianTermsGatherPass1Kernel(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        double* __restrict__ scratch,
        const ProcessSetParam devstruct)
    {
      const std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;
      const std::size_t VL = state_variables_param.vector_length_;
      const std::size_t number_of_groups = (number_of_grid_cells + VL - 1) / VL;
      const std::size_t group_id = tid / VL;
      const std::size_t local_tid = tid % VL;

      if (tid >= number_of_grid_cells)
        return;

      const double* __restrict__ d_rate_constants =
          rate_constants_param.d_data_ + group_id * (rate_constants_param.number_of_elements_ / number_of_groups);
      const double* __restrict__ d_state_variables =
          state_variables_param.d_data_ + group_id * (state_variables_param.number_of_elements_ / number_of_groups);
      double* __restrict__ my_scratch = scratch + group_id * number_of_process_infos * VL;

      std::size_t react_off = 0;
      for (std::size_t i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
      {
        const ProcessInfoParam pi = devstruct.jacobian_process_info_[i_proc];
        double d_rate = d_rate_constants[pi.process_id_ * VL + local_tid];
        for (std::size_t r = 0; r < pi.number_of_dependent_reactants_; ++r)
          d_rate *= d_state_variables[devstruct.jacobian_reactant_ids_[react_off + r] * VL + local_tid];
        my_scratch[i_proc * VL + local_tid] = d_rate;
        react_off += pi.number_of_dependent_reactants_;
      }
    }  // end of SubtractJacobianTermsGatherPass1Kernel

    /// Two-pass gather, Pass 2: gather from scratch using the CSR and write each Jacobian
    /// entry exactly once. Scratch layout identical to Pass 1.
    __global__ void SubtractJacobianTermsGatherPass2Kernel(
        const double* __restrict__ scratch,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      const std::size_t cell_id = blockIdx.x * blockDim.x + threadIdx.x;
      const std::size_t jac_idx = blockIdx.y * blockDim.y + threadIdx.y;

      const std::size_t number_of_grid_cells = jacobian_param.number_of_grid_cells_;
      if (cell_id >= number_of_grid_cells || jac_idx >= devstruct.jac_gather_entries_size_)
        return;

      const std::size_t VL = state_variables_param.vector_length_;
      const std::size_t number_of_groups = (number_of_grid_cells + VL - 1) / VL;
      const std::size_t group_id = cell_id / VL;
      const std::size_t local_tid = cell_id % VL;
      const std::size_t N_proc = devstruct.jacobian_process_info_size_;

      const double* __restrict__ my_scratch = scratch + group_id * N_proc * VL;
      double* __restrict__ d_jacobian =
          jacobian_param.d_data_ + group_id * (jacobian_param.number_of_elements_ / number_of_groups);

      const std::size_t start = devstruct.jac_gather_offsets_[jac_idx];
      const std::size_t end = devstruct.jac_gather_offsets_[jac_idx + 1];

      double total = 0.0;
      for (std::size_t s = start; s < end; ++s)
        total += devstruct.jac_gather_coeffs_[s] * my_scratch[devstruct.jac_gather_proc_idx_[s] * VL + local_tid];

      // Single write per (jac_idx, local_tid) — no atomic needed
      d_jacobian[devstruct.jac_gather_unique_flat_ids_[jac_idx] + local_tid] += total;
    }  // end of SubtractJacobianTermsGatherPass2Kernel

    void SubtractJacobianTermsTwoPassGatherDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      const std::size_t VL = state_variables_param.vector_length_;
      const std::size_t N_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t N_groups = (N_cells + VL - 1) / VL;
      const std::size_t N_proc = devstruct.jacobian_process_info_size_;

      // Scratch: one d_rate_d_ind per (group, process, lane)
      double* scratch = nullptr;
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&scratch, N_groups * N_proc * VL * sizeof(double), stream), "cudaMalloc scratch");

      // Pass 1: compute and store d_rate_d_ind for every process
      const std::size_t nblocks1 = (N_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SubtractJacobianTermsGatherPass1Kernel<<<nblocks1, BLOCK_SIZE, 0, stream>>>(
          rate_constants_param, state_variables_param, scratch, devstruct);

      // Pass 2: gather weighted contributions to each Jacobian entry
      constexpr std::size_t BLOCK_DIM_X = 32;
      constexpr std::size_t BLOCK_DIM_Y = 4;
      static_assert(BLOCK_DIM_X * BLOCK_DIM_Y == BLOCK_SIZE, "2D block must equal BLOCK_SIZE");
      const std::size_t N_unique = devstruct.jac_gather_entries_size_;
      const dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
      const dim3 grid((N_cells + BLOCK_DIM_X - 1) / BLOCK_DIM_X, (N_unique + BLOCK_DIM_Y - 1) / BLOCK_DIM_Y);
      SubtractJacobianTermsGatherPass2Kernel<<<grid, block, 0, stream>>>(
          scratch, state_variables_param, jacobian_param, devstruct);

      CHECK_CUDA_ERROR(cudaFreeAsync(scratch, stream), "cudaFree scratch");
    }  // end of SubtractJacobianTermsTwoPassGatherDriver

  }  // namespace cuda
}  // namespace micm
