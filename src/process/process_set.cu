// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that calculates the forcing terms on the device
    __global__ void AddForcingTermsKernel(
        double* d_rate_constants,
        double* d_state_variables,
        double* d_forcing,
        ProcessSetParam devstruct,
        size_t n_grids,
        size_t n_reactions,
        size_t n_species)
    {
      /// Local device variables
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      size_t react_id_offset, prod_id_offset, yield_offset;
      size_t* d_number_of_reactants = devstruct.number_of_reactants_;
      size_t* d_reactant_ids = devstruct.reactant_ids_;
      size_t* d_number_of_products = devstruct.number_of_products_;
      size_t* d_product_ids = devstruct.product_ids_;
      double* d_yields = devstruct.yields_;

      if (tid < n_grids)
      {
        react_id_offset = 0;
        prod_id_offset = 0;
        yield_offset = 0;
        for (std::size_t i_rxn = 0; i_rxn < n_reactions; ++i_rxn)
        {
          double rate = d_rate_constants[i_rxn * n_grids + tid];
          for (std::size_t i_react = 0; i_react < d_number_of_reactants[i_rxn]; ++i_react)
            rate *= d_state_variables[d_reactant_ids[react_id_offset + i_react] * n_grids + tid];
          for (std::size_t i_react = 0; i_react < d_number_of_reactants[i_rxn]; ++i_react)
          {
            d_forcing[d_reactant_ids[react_id_offset + i_react] * n_grids + tid] -= rate;
          }
          for (std::size_t i_prod = 0; i_prod < d_number_of_products[i_rxn]; ++i_prod)
          {
            size_t index = d_product_ids[prod_id_offset + i_prod] * n_grids + tid;
            d_forcing[index] += d_yields[yield_offset + i_prod] * rate;
          }
          react_id_offset += d_number_of_reactants[i_rxn];
          prod_id_offset += d_number_of_products[i_rxn];
          yield_offset += d_number_of_products[i_rxn];
        }  // end of loop over number of reactions
      }    // end of checking a valid CUDA thread id
    }      // end of AddForcingTerms_kernel

    /// This is the CUDA kernel that forms the Jacobian matrix on the device
    __global__ void AddJacobianTermsKernel(
        double* d_rate_constants,
        double* d_state_variables,
        double* d_jacobian,
        ProcessSetParam devstruct,
        size_t n_grids,
        size_t n_reactions)
    {
      /// Local device variables
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      size_t react_ids_offset = 0;
      size_t yields_offset = 0;
      size_t flat_id_offset = 0;
      size_t* d_number_of_reactants = devstruct.number_of_reactants_;
      size_t* d_reactant_ids = devstruct.reactant_ids_;
      size_t* d_number_of_products = devstruct.number_of_products_;
      size_t* d_jacobian_flat_ids = devstruct.jacobian_flat_ids_;
      double* d_yields = devstruct.yields_;

      if (tid < n_grids)
      {
        // loop over reactions in a grid
        for (size_t i_rxn = 0; i_rxn < n_reactions; ++i_rxn)
        {
          // loop over reactants in a reaction
          for (size_t i_ind = 0; i_ind < d_number_of_reactants[i_rxn]; ++i_ind)
          {
            double d_rate_d_ind = d_rate_constants[i_rxn * n_grids + tid];
            for (size_t i_react = 0; i_react < d_number_of_reactants[i_rxn]; ++i_react)
            {
              if (i_react != i_ind)
              {
                d_rate_d_ind *= d_state_variables[d_reactant_ids[react_ids_offset + i_react] * n_grids + tid];
              }
            }
            for (size_t i_dep = 0; i_dep < d_number_of_reactants[i_rxn]; ++i_dep)
            {
              size_t jacobian_idx = d_jacobian_flat_ids[flat_id_offset] + tid;
              d_jacobian[jacobian_idx] -= d_rate_d_ind;
              flat_id_offset++;
            }
            for (size_t i_dep = 0; i_dep < d_number_of_products[i_rxn]; ++i_dep)
            {
              size_t jacobian_idx = d_jacobian_flat_ids[flat_id_offset] + tid;
              d_jacobian[jacobian_idx] += d_yields[yields_offset + i_dep] * d_rate_d_ind;
              flat_id_offset++;
            }
          }  // loop over reactants in a reaction
          react_ids_offset += d_number_of_reactants[i_rxn];
          yields_offset += d_number_of_products[i_rxn];
        }  // end of loop over reactions in a grid
      }    // end of checking a CUDA thread id
    }      // end of AddJacobianTermsKernel

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

      /// Allocate memory space on the device
      cudaMalloc(&(devstruct.number_of_reactants_), number_of_reactants_bytes);
      cudaMalloc(&(devstruct.reactant_ids_), reactant_ids_bytes);
      cudaMalloc(&(devstruct.number_of_products_), number_of_products_bytes);
      cudaMalloc(&(devstruct.product_ids_), product_ids_bytes);
      cudaMalloc(&(devstruct.yields_), yields_bytes);

      /// Copy the data from host to device
      cudaMemcpy(
          devstruct.number_of_reactants_,
          hoststruct.number_of_reactants_,
          number_of_reactants_bytes,
          cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.reactant_ids_, hoststruct.reactant_ids_, reactant_ids_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(
          devstruct.number_of_products_, hoststruct.number_of_products_, number_of_products_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.product_ids_, hoststruct.product_ids_, product_ids_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.yields_, hoststruct.yields_, yields_bytes, cudaMemcpyHostToDevice);

      devstruct.number_of_reactants_size_ = hoststruct.number_of_reactants_size_;
      devstruct.reactant_ids_size_ = hoststruct.reactant_ids_size_;
      devstruct.number_of_products_size_ = hoststruct.number_of_products_size_;
      devstruct.product_ids_size_ = hoststruct.product_ids_size_;
      devstruct.yields_size_ = hoststruct.yields_size_;

      return devstruct;
    }

    /// This is the function that will copy the jacobian_flat_id
    ///   of class "ProcessSet" to the device, after the matrix
    ///   structure is known;
    void CopyJacobiFlatId(ProcessSetParam& hoststruct, ProcessSetParam& devstruct)
    {
      /// Calculate the memory space
      size_t jacobian_flat_ids_bytes = sizeof(size_t) * hoststruct.jacobian_flat_ids_size_;

      /// Allocate memory space on the device
      cudaMalloc(&(devstruct.jacobian_flat_ids_), jacobian_flat_ids_bytes);

      /// Copy the data from host to device
      cudaMemcpy(
          devstruct.jacobian_flat_ids_, hoststruct.jacobian_flat_ids_, jacobian_flat_ids_bytes, cudaMemcpyHostToDevice);

      devstruct.jacobian_flat_ids_size_ = hoststruct.jacobian_flat_ids_size_;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaProcessSet" on the device
    void FreeConstData(ProcessSetParam& devstruct)
    {
      cudaFree(devstruct.number_of_reactants_);
      cudaFree(devstruct.reactant_ids_);
      cudaFree(devstruct.number_of_products_);
      cudaFree(devstruct.product_ids_);
      cudaFree(devstruct.yields_);
      if (devstruct.jacobian_flat_ids_ != nullptr)
        cudaFree(devstruct.jacobian_flat_ids_);
    }

    std::chrono::nanoseconds AddJacobianTermsKernelDriver(
        CudaMatrixParam& matrixParam,
        CudaSparseMatrixParam& sparseMatrix,
        const ProcessSetParam& devstruct)
    {
      // create device pointers
      double* d_rate_constants;
      double* d_state_variables;
      double* d_jacobian;

      // allocate device memory
      cudaMalloc(&d_rate_constants, sizeof(double) * matrixParam.n_grids_ * matrixParam.n_reactions_);
      cudaMalloc(&d_state_variables, sizeof(double) * matrixParam.n_grids_ * matrixParam.n_species_);
      cudaMalloc(&d_jacobian, sizeof(double) * sparseMatrix.jacobian_size_);

      // transfer data from host to device
      cudaMemcpy(
          d_rate_constants,
          matrixParam.rate_constants_,
          sizeof(double) * matrixParam.n_grids_ * matrixParam.n_reactions_,
          cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_state_variables,
          matrixParam.state_variables_,
          sizeof(double) * matrixParam.n_grids_ * matrixParam.n_species_,
          cudaMemcpyHostToDevice);
      cudaMemcpy(d_jacobian, sparseMatrix.jacobian_, sizeof(double) * sparseMatrix.jacobian_size_, cudaMemcpyHostToDevice);

      // setup kernel
      size_t num_blocks = (matrixParam.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      size_t n_reactions = matrixParam.n_reactions_;
      size_t n_grids = matrixParam.n_grids_;

      // launch kernel and measure time performance
      auto startTime = std::chrono::high_resolution_clock::now();
      AddJacobianTermsKernel<<<num_blocks, BLOCK_SIZE>>>(
          d_rate_constants, d_state_variables, d_jacobian, devstruct, n_grids, n_reactions);
      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

      // copy the data from device to host
      cudaMemcpy(sparseMatrix.jacobian_, d_jacobian, sizeof(double) * sparseMatrix.jacobian_size_, cudaMemcpyDeviceToHost);

      // clean up
      cudaFree(d_rate_constants);
      cudaFree(d_state_variables);
      cudaFree(d_jacobian);

      return kernel_duration;
    }  // end of AddJacobian_kernelSetup

    std::chrono::nanoseconds AddForcingTermsKernelDriver(CudaMatrixParam& matrixParam, const ProcessSetParam& devstruct)
    {
      // device pointer to vectorss
      double* d_rate_constants;
      double* d_state_variables;
      double* d_forcing;

      // allocate device memory
      cudaMalloc(&d_rate_constants, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_reactions_));
      cudaMalloc(&d_state_variables, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_));
      cudaMalloc(&d_forcing, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_));

      // copy data from host memory to device memory
      cudaMemcpy(
          d_rate_constants,
          matrixParam.rate_constants_,
          sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_reactions_),
          cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_state_variables,
          matrixParam.state_variables_,
          sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_),
          cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_forcing,
          matrixParam.forcing_,
          sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_),
          cudaMemcpyHostToDevice);

      int num_block = (matrixParam.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      size_t n_grids = matrixParam.n_grids_;
      size_t n_reactions = matrixParam.n_reactions_;
      size_t n_species = matrixParam.n_species_;

      // launch kernel and measure time performance
      auto startTime = std::chrono::high_resolution_clock::now();
      AddForcingTermsKernel<<<num_block, BLOCK_SIZE>>>(
          d_rate_constants, d_state_variables, d_forcing, devstruct, n_grids, n_reactions, n_species);
      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

      // copy data from device memory to host memory
      cudaMemcpy(matrixParam.forcing_, d_forcing, sizeof(double) * (n_grids * n_species), cudaMemcpyDeviceToHost);

      // clean up
      cudaFree(d_rate_constants);
      cudaFree(d_state_variables);
      cudaFree(d_forcing);

      return kernel_duration;
    }  // end of AddForcingTerms_kernelSetup
  }    // namespace cuda
}  // namespace micm
