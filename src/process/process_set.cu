// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <micm/util/cuda_param.hpp>

// device pointers passing to AddJacobianTermsKernel()
typedef struct JacobianDevice
{
  double* rate_constants_;
  double* state_variables_;
  double* jacobian_;
  size_t* number_of_reactants_;
  size_t* reactant_ids_;
  size_t* number_of_products_;
  double* yields_;
  size_t* jacobian_flat_ids_;
};

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that calculates the forcing terms on the device
    __global__ void AddForcingTermsKernel(double* d_rate_constants, double* d_state_variables,
                                          double* d_forcing, ProcessSetParam devstruct,
                                          size_t n_grids, size_t n_reactions, size_t n_species)
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
        }  // for loop over number of reactions
      }    // if check for valid Cuda threads
    }      // end of AddForcingTerms_kernel

    __global__ void FormJacobianMatrixKernel(JacobianDevice* device, size_t n_grids, size_t n_reactions)
    {
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      size_t react_ids_offset = 0;
      size_t yields_offset = 0;
      size_t flat_id_offset = 0;
      size_t* number_of_reactants = device->number_of_reactants_;
      size_t* jacobian_flat_ids = device->jacobian_flat_ids_;
      size_t* number_of_products = device->number_of_products_;
      double* jacobian = device->jacobian_;

      if (tid < n_grids)
      {
        // loop over reactions in a grid
        for (size_t i_rxn = 0; i_rxn < n_reactions; ++i_rxn)
        {
          // loop over reactants in a reaction
          for (size_t i_ind = 0; i_ind < number_of_reactants[i_rxn]; ++i_ind)
          {
            double d_rate_d_ind = device->rate_constants_[i_rxn * n_grids + tid];
            for (size_t i_react = 0; i_react < number_of_reactants[i_rxn]; ++i_react)
            {
              if (i_react != i_ind)
              {
                d_rate_d_ind *= device->state_variables_[device->reactant_ids_[react_ids_offset + i_react] * n_grids + tid];
              }
            }
            for (size_t i_dep = 0; i_dep < number_of_reactants[i_rxn]; ++i_dep)
            {
              size_t jacobian_idx = jacobian_flat_ids[flat_id_offset] + tid;
              jacobian[jacobian_idx] -= d_rate_d_ind;
              flat_id_offset++;
            }
            for (size_t i_dep = 0; i_dep < number_of_products[i_rxn]; ++i_dep)
            {
              size_t jacobian_idx = jacobian_flat_ids[flat_id_offset] + tid;
              jacobian[jacobian_idx] += device->yields_[yields_offset + i_dep] * d_rate_d_ind;
              flat_id_offset++;
            }
          }  // loop over reactants in a reaction
          react_ids_offset += number_of_reactants[i_rxn];
          yields_offset += number_of_products[i_rxn];
        }  // loop over reactions in a grid
      }    // check valid tid
    }      // end of AddJacobianTerms_kernel

    /// This is the function that will copy the constant data
    ///   members of class "CudaProcessSet" to the device
    ProcessSetParam CopyConstData(ProcessSetParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      size_t number_of_reactants_bytes = sizeof(size_t) * hoststruct.number_of_reactants_size_;
      size_t reactant_ids_bytes = sizeof(size_t) * hoststruct.reactant_ids_size_;
      size_t number_of_products_bytes = sizeof(size_t) * hoststruct.number_of_products_size_;
      size_t product_ids_bytes = sizeof(size_t) * hoststruct.product_ids_size_;
      size_t yields_bytes = sizeof(double) * hoststruct.yields_size_;
      size_t jacobian_flat_ids_bytes = sizeof(size_t) * hoststruct.jacobian_flat_ids_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      ProcessSetParam devstruct;

      /// Allocate memory space on the device
      cudaMalloc(&(devstruct.number_of_reactants_), number_of_reactants_bytes);
      cudaMalloc(&(devstruct.reactant_ids_), reactant_ids_bytes);
      cudaMalloc(&(devstruct.number_of_products_), number_of_products_bytes);
      cudaMalloc(&(devstruct.product_ids_), product_ids_bytes);
      cudaMalloc(&(devstruct.yields_), yields_bytes);
      cudaMalloc(&(devstruct.jacobian_flat_ids_), jacobian_flat_ids_bytes);

      /// Copy the data from host to device
      cudaMemcpy(devstruct.number_of_reactants_, hoststruct.number_of_reactants_, number_of_reactants_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.reactant_ids_, hoststruct.reactant_ids_, reactant_ids_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.number_of_products_, hoststruct.number_of_products_, number_of_products_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.product_ids_, hoststruct.product_ids_, product_ids_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.yields_, hoststruct.yields_, yields_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.jacobian_flat_ids_, hoststruct.jacobian_flat_ids_, jacobian_flat_ids_bytes, cudaMemcpyHostToDevice);

      devstruct.number_of_reactants_size_ = hoststruct.number_of_reactants_size_;
      devstruct.reactant_ids_size_ = hoststruct.reactant_ids_size_;
      devstruct.number_of_products_size_ = hoststruct.number_of_products_size_;
      devstruct.product_ids_size_ = hoststruct.product_ids_size_;
      devstruct.yields_size_ = hoststruct.yields_size_;
      devstruct.jacobian_flat_ids_size_ = hoststruct.jacobian_flat_ids_size_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaProcessSet" on the device
    void FreeConstData(ProcessSetParam& devstruct)
    {
      cudaFree(&(devstruct.number_of_reactants_));
      cudaFree(&(devstruct.reactant_ids_));
      cudaFree(&(devstruct.number_of_products_));
      cudaFree(&(devstruct.product_ids_));
      cudaFree(&(devstruct.yields_));
      cudaFree(&(devstruct.jacobian_flat_ids_));
    }

    std::chrono::nanoseconds FormJacobianMatrixKernelDriver(
        CudaMatrixParam& matrixParam,
        CudaSparseMatrixParam& sparseMatrix,
        ProcessSetParam& processSet)
    {
      // create device pointers
      double* d_rate_constants;
      double* d_state_variables;
      double* d_jacobian;
      size_t* d_number_of_reactants;
      size_t* d_reactant_ids;
      size_t* d_number_of_products;
      double* d_yields;
      size_t* d_jacobian_flat_ids;
      JacobianDevice* device;

      // allocate device memory
      cudaMalloc(&d_rate_constants, sizeof(double) * matrixParam.n_grids_ * matrixParam.n_reactions_);
      cudaMalloc(&d_state_variables, sizeof(double) * matrixParam.n_grids_ * matrixParam.n_species_);
      cudaMalloc(&d_jacobian, sizeof(double) * sparseMatrix.jacobian_size_);
      cudaMalloc(&d_number_of_reactants, sizeof(size_t) * matrixParam.n_reactions_);
      cudaMalloc(&d_reactant_ids, sizeof(size_t) * processSet.reactant_ids_size_);
      cudaMalloc(&d_number_of_products, sizeof(size_t) * matrixParam.n_reactions_);
      cudaMalloc(&d_yields, sizeof(double) * processSet.yields_size_);
      cudaMalloc(&d_jacobian_flat_ids, sizeof(size_t) * processSet.jacobian_flat_ids_size_);
      cudaMalloc(&device, sizeof(JacobianDevice));

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
      cudaMemcpy(
          d_number_of_reactants,
          processSet.number_of_reactants_,
          sizeof(size_t) * matrixParam.n_reactions_,
          cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_reactant_ids, processSet.reactant_ids_, sizeof(size_t) * processSet.reactant_ids_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_number_of_products,
          processSet.number_of_products_,
          sizeof(size_t) * matrixParam.n_reactions_,
          cudaMemcpyHostToDevice);
      cudaMemcpy(d_yields, processSet.yields_, sizeof(double) * processSet.yields_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_jacobian_flat_ids,
          processSet.jacobian_flat_ids_,
          sizeof(size_t) * processSet.jacobian_flat_ids_size_,
          cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->rate_constants_), &d_rate_constants, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->state_variables_), &d_state_variables, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->jacobian_), &d_jacobian, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->number_of_reactants_), &d_number_of_reactants, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->reactant_ids_), &d_reactant_ids, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->number_of_products_), &d_number_of_products, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->yields_), &d_yields, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->jacobian_flat_ids_), &d_jacobian_flat_ids, sizeof(size_t*), cudaMemcpyHostToDevice);

      // setup kernel
      size_t total_blocks = (matrixParam.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;

      size_t n_reactions = matrixParam.n_reactions_;
      size_t n_grids = matrixParam.n_grids_;
      // launch kernel and measure time performance
      auto startTime = std::chrono::high_resolution_clock::now();
      FormJacobianMatrixKernel<<<total_blocks, BLOCK_SIZE>>>(device, n_grids, n_reactions);
      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

      cudaMemcpy(sparseMatrix.jacobian_, d_jacobian, sizeof(double) * sparseMatrix.jacobian_size_, cudaMemcpyDeviceToHost);
      // clean up
      cudaFree(d_rate_constants);
      cudaFree(d_state_variables);
      cudaFree(d_jacobian);
      cudaFree(d_number_of_reactants);
      cudaFree(d_reactant_ids);
      cudaFree(d_number_of_products);
      cudaFree(d_yields);
      cudaFree(d_jacobian_flat_ids);
      cudaFree(device);
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
      cudaMemcpy(d_rate_constants, matrixParam.rate_constants_, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_reactions_), cudaMemcpyHostToDevice);
      cudaMemcpy(d_state_variables, matrixParam.state_variables_, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_), cudaMemcpyHostToDevice);
      cudaMemcpy(d_forcing, matrixParam.forcing_, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_), cudaMemcpyHostToDevice);

      int num_block = (matrixParam.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      size_t n_grids = matrixParam.n_grids_;
      size_t n_reactions = matrixParam.n_reactions_;
      size_t n_species = matrixParam.n_species_;

      // launch kernel and measure time performance
      auto startTime = std::chrono::high_resolution_clock::now();
      AddForcingTermsKernel<<<num_block, BLOCK_SIZE>>>(d_rate_constants, d_state_variables, d_forcing,
                                                       devstruct, n_grids, n_reactions, n_species);
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