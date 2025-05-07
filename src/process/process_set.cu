// Copyright (C) 2023-2025 University Corporation for Atmospheric Research-National Center for Atmospheric Research
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
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      size_t react_id_offset, prod_id_offset, yield_offset;
      const size_t* const __restrict__ d_number_of_reactants = devstruct.number_of_reactants_;
      const size_t* const __restrict__ d_reactant_ids = devstruct.reactant_ids_;
      const size_t* const __restrict__ d_number_of_products = devstruct.number_of_products_;
      const size_t* const __restrict__ d_product_ids = devstruct.product_ids_;
      const double* const __restrict__ d_yields = devstruct.yields_;
      const size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const size_t number_of_reactions =
          rate_constants_param.number_of_elements_ / rate_constants_param.number_of_grid_cells_;
      const double* const __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* const __restrict__ d_state_variables = state_variables_param.d_data_;
      double* const __restrict__ d_forcing = forcing_param.d_data_;

      if (tid < number_of_grid_cells)
      {
        react_id_offset = 0;
        prod_id_offset = 0;
        yield_offset = 0;
        for (std::size_t i_rxn = 0; i_rxn < number_of_reactions; ++i_rxn)
        {
          double rate = d_rate_constants[(i_rxn * number_of_grid_cells) + tid];
          for (std::size_t i_react = 0; i_react < d_number_of_reactants[i_rxn]; ++i_react)
          {
            rate *= d_state_variables[(d_reactant_ids[react_id_offset + i_react] * number_of_grid_cells) + tid];
          }
          for (std::size_t i_react = 0; i_react < d_number_of_reactants[i_rxn]; ++i_react)
          {
            d_forcing[(d_reactant_ids[react_id_offset + i_react] * number_of_grid_cells) + tid] -= rate;
          }
          for (std::size_t i_prod = 0; i_prod < d_number_of_products[i_rxn]; ++i_prod)
          {
            size_t index = d_product_ids[prod_id_offset + i_prod] * number_of_grid_cells + tid;
            d_forcing[index] += d_yields[yield_offset + i_prod] * rate;
          }
          react_id_offset += d_number_of_reactants[i_rxn];
          prod_id_offset += d_number_of_products[i_rxn];
          yield_offset += d_number_of_products[i_rxn];
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
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      /// Local device variables
      size_t react_ids_offset = 0;
      size_t yields_offset = 0;
      size_t flat_id_offset = 0;
      const ProcessInfoParam* const __restrict__ d_jacobian_process_info = devstruct.jacobian_process_info_;
      const size_t* const __restrict__ d_reactant_ids = devstruct.jacobian_reactant_ids_;
      const double* const __restrict__ d_yields = devstruct.jacobian_yields_;
      const size_t* const __restrict__ d_jacobian_flat_ids = devstruct.jacobian_flat_ids_;
      const size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const size_t number_of_process_infos = devstruct.jacobian_process_info_size_;
      const double* const __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* const __restrict__ d_state_variables = state_variables_param.d_data_;
      double* const __restrict__ d_jacobian = jacobian_param.d_data_;

      if (tid < number_of_grid_cells)
      {
        // loop over reactions in a grid
        for (size_t i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
        {
          const ProcessInfoParam& process_info = d_jacobian_process_info[i_proc];
          // Calculate d_rate/d_ind
          double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * number_of_grid_cells) + tid];
          for (size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
          {
            d_rate_d_ind *= d_state_variables[(d_reactant_ids[react_ids_offset + i_react] * number_of_grid_cells) + tid];
          }
          for (size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
          {
            size_t jacobian_idx = d_jacobian_flat_ids[flat_id_offset] + tid;
            d_jacobian[jacobian_idx] += d_rate_d_ind;
            flat_id_offset++;
          }
          for (size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
          {
            size_t jacobian_idx = d_jacobian_flat_ids[flat_id_offset] + tid;
            d_jacobian[jacobian_idx] -= d_yields[yields_offset + i_dep] * d_rate_d_ind;
            flat_id_offset++;
          }
          react_ids_offset += process_info.number_of_dependent_reactants_;
          yields_offset += process_info.number_of_products_;
        }  // end of loop over reactions in a grid cell
      }  // end of checking a CUDA thread id
    }  // end of SubtractJacobianTermsKernel

    /// This is the function that will copy the constant data
    ///   members of class "CudaProcessSet" to the device,
    ///   except for the "jacobian_flat_id" because it is unknown now
    ProcessSetParam CopyConstData(ProcessSetParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      size_t number_of_reactants_bytes = sizeof(size_t) * hoststruct.number_of_reactants_size_;
      size_t reactant_ids_bytes = sizeof(size_t) * hoststruct.reactant_ids_size_;
      size_t number_of_products_bytes = sizeof(size_t) * hoststruct.number_of_products_size_;
      size_t product_ids_bytes = sizeof(size_t) * hoststruct.product_ids_size_;
      size_t yields_bytes = sizeof(double) * hoststruct.yields_size_;

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
      size_t jacobian_process_info_bytes = sizeof(ProcessInfoParam) * hoststruct.jacobian_process_info_size_;
      size_t jacobian_reactant_ids_bytes = sizeof(size_t) * hoststruct.jacobian_reactant_ids_size_;
      size_t jacobian_product_ids_bytes = sizeof(size_t) * hoststruct.jacobian_product_ids_size_;
      size_t jacobian_yields_bytes = sizeof(double) * hoststruct.jacobian_yields_size_;
      size_t jacobian_flat_ids_bytes = sizeof(size_t) * hoststruct.jacobian_flat_ids_size_;

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
    }

    void SubtractJacobianTermsKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
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
      size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      AddForcingTermsKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, forcing_param, devstruct);
    }  // end of AddForcingTermsKernelDriver
  }  // namespace cuda
}  // namespace micm
