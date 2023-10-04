// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <iostream>
#include <micm/util/cuda_param.hpp>

// device pointers passing to AddForcingTermsKernel()
typedef struct forcingDevice
{
  double* rate_constants;
  double* state_variables;
  double* forcing;
  size_t* number_of_reactants;
  size_t* reactant_ids;
  size_t* number_of_products;
  size_t* product_ids;
  double* yields;
};
// device pointers passing to AddJacobianTermsKernel()
typedef struct jacobianDevice
{
  double* rate_constants;
  double* state_variables;
  double* jacobian;
  size_t* number_of_reactants;
  size_t* reactant_ids;
  size_t* number_of_products;
  double* yields;
  size_t* jacobian_flat_ids;
};

namespace micm
{
  namespace cuda
  {
    // flipped memory layout
    __global__ void AddForcingTermsKernel(forcingDevice* device, size_t n_grids, size_t n_reactions, size_t n_species)

    {
      // define thread index
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      size_t react_id_offset, prod_id_offset, yield_offset;
      double* forcing = device->forcing;
      size_t* number_of_reactants = device->number_of_reactants;
      size_t* reactant_ids = device->reactant_ids;
      size_t* number_of_products = device->number_of_products;
      if (tid < n_grids)
      {
        react_id_offset = 0;
        prod_id_offset = 0;
        yield_offset = 0;
        for (std::size_t i_rxn = 0; i_rxn < n_reactions; ++i_rxn)
        {
          double rate = device->rate_constants[i_rxn * n_grids + tid];
          for (std::size_t i_react = 0; i_react < number_of_reactants[i_rxn]; ++i_react)
            rate *= device->state_variables[reactant_ids[react_id_offset + i_react] * n_grids + tid];
          for (std::size_t i_react = 0; i_react < number_of_reactants[i_rxn]; ++i_react)
          {
            forcing[reactant_ids[react_id_offset + i_react] * n_grids + tid] -= rate;
          }
          for (std::size_t i_prod = 0; i_prod < number_of_products[i_rxn]; ++i_prod)
          {
            size_t index = device->product_ids[prod_id_offset + i_prod] * n_grids + tid;
            forcing[index] += device->yields[yield_offset + i_prod] * rate;
          }
          react_id_offset += number_of_reactants[i_rxn];
          prod_id_offset += number_of_products[i_rxn];
          yield_offset += number_of_products[i_rxn];
        }  // for loop over number of reactions
      }    // if check for valid Cuda threads
    }      // end of AddForcingTerms_kernel

    __global__ void AddJacobianTermsKernel(jacobianDevice* device, size_t n_grids, size_t n_reactions)
    {
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      size_t react_ids_offset = 0;
      size_t yields_offset = 0;
      size_t flat_id_offset = 0;
      size_t* number_of_reactants = device->number_of_reactants;
      size_t* jacobian_flat_ids = device->jacobian_flat_ids;
      size_t* number_of_products = device->number_of_products;
      double* jacobian = device->jacobian;

      if (tid < n_grids)
      {
        // loop over reactions in a grid
        for (size_t i_rxn = 0; i_rxn < n_reactions; ++i_rxn)
        {
          // loop over reactants in a reaction
          for (size_t i_ind = 0; i_ind < number_of_reactants[i_rxn]; ++i_ind)
          {
            double d_rate_d_ind = device->rate_constants[i_rxn * n_grids + tid];
            for (size_t i_react = 0; i_react < number_of_reactants[i_rxn]; ++i_react)
            {
              if (i_react != i_ind)
              {
                d_rate_d_ind *= device->state_variables[device->reactant_ids[react_ids_offset + i_react] * n_grids + tid];
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
              jacobian[jacobian_idx] += device->yields[yields_offset + i_dep] * d_rate_d_ind;
              flat_id_offset++;
            }
          }  // loop over reactants in a reaction
          react_ids_offset += number_of_reactants[i_rxn];
          yields_offset += number_of_products[i_rxn];
        }  // loop over reactions in a grid
      }    // check valid tid
    }      // end of AddJacobianTerms_kernel

    std::chrono::nanoseconds AddJacobianTermsKernelDriver(
        CudaMatrixParam& matrixParam,
        CudaSparseMatrixParam& sparseMatrix,
        CudaProcessSetParam& processSet)
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
      jacobianDevice* device;

      // allocate device memory
      cudaMalloc(&d_rate_constants, sizeof(double) * matrixParam.n_grids_ * matrixParam.n_reactions_);
      cudaMalloc(&d_state_variables, sizeof(double) * matrixParam.n_grids_ * matrixParam.n_species_);
      cudaMalloc(&d_jacobian, sizeof(double) * sparseMatrix.jacobian_size_);
      cudaMalloc(&d_number_of_reactants, sizeof(size_t) * matrixParam.n_reactions_);
      cudaMalloc(&d_reactant_ids, sizeof(size_t) * processSet.reactant_ids_size_);
      cudaMalloc(&d_number_of_products, sizeof(size_t) * matrixParam.n_reactions_);
      cudaMalloc(&d_yields, sizeof(double) * processSet.yields_size_);
      cudaMalloc(&d_jacobian_flat_ids, sizeof(size_t) * processSet.jacobian_flat_ids_size_);
      cudaMalloc(&device, sizeof(jacobianDevice));

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
      cudaMemcpy(&(device->rate_constants), &d_rate_constants, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->state_variables), &d_state_variables, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->jacobian), &d_jacobian, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->number_of_reactants), &d_number_of_reactants, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->reactant_ids), &d_reactant_ids, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->number_of_products), &d_number_of_products, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->yields), &d_yields, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->jacobian_flat_ids), &d_jacobian_flat_ids, sizeof(size_t*), cudaMemcpyHostToDevice);

      // setup kernel
      size_t total_blocks = (matrixParam.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;

      size_t n_reactions = matrixParam.n_reactions_;
      size_t n_grids = matrixParam.n_grids_;
      // launch kernel and measure time performance
      auto startTime = std::chrono::high_resolution_clock::now();
      AddJacobianTermsKernel<<<total_blocks, BLOCK_SIZE>>>(device, n_grids, n_reactions);
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

    std::chrono::nanoseconds AddForcingTermsKernelDriver(CudaMatrixParam& matrixParam, CudaProcessSetParam& processSet)
    {
      // device pointer to vectorss
      double* d_rate_constants;
      double* d_state_variables;
      double* d_forcing;
      double* d_yields;
      size_t* d_number_of_reactants;
      size_t* d_reactant_ids;
      size_t* d_number_of_products;
      size_t* d_product_ids;
      forcingDevice* device;

      // allocate device memory
      cudaMalloc(&d_rate_constants, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_reactions_));
      cudaMalloc(&d_state_variables, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_));
      cudaMalloc(&d_forcing, sizeof(double) * (matrixParam.n_grids_ * matrixParam.n_species_));
      cudaMalloc(&d_number_of_reactants, sizeof(size_t) * matrixParam.n_reactions_);
      cudaMalloc(&d_reactant_ids, sizeof(size_t) * processSet.reactant_ids_size_);
      cudaMalloc(&d_number_of_products, sizeof(size_t) * matrixParam.n_reactions_);
      cudaMalloc(&d_product_ids, sizeof(size_t) * processSet.product_ids_size_);
      cudaMalloc(&d_yields, sizeof(double) * processSet.yields_size_);
      cudaMalloc(&device, sizeof(forcingDevice));

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
      cudaMemcpy(
          d_product_ids, processSet.product_ids_, sizeof(size_t) * processSet.product_ids_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(d_yields, processSet.yields_, sizeof(double) * processSet.yields_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->rate_constants), &d_rate_constants, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->state_variables), &d_state_variables, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->forcing), &d_forcing, sizeof(double*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->number_of_reactants), &d_number_of_reactants, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->reactant_ids), &d_reactant_ids, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->number_of_products), &d_number_of_products, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->product_ids), &d_product_ids, sizeof(size_t*), cudaMemcpyHostToDevice);
      cudaMemcpy(&(device->yields), &d_yields, sizeof(double*), cudaMemcpyHostToDevice);

      // total thread count == number of grid cells

      int num_block = (matrixParam.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;

      size_t n_grids = matrixParam.n_grids_;
      size_t n_reactions = matrixParam.n_reactions_;
      size_t n_species = matrixParam.n_species_;

      // launch kernel and measure time performance
      auto startTime = std::chrono::high_resolution_clock::now();
      AddForcingTermsKernel<<<num_block, BLOCK_SIZE>>>(device, n_grids, n_reactions, n_species);
      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

      // copy data from device memory to host memory
      cudaMemcpy(matrixParam.forcing_, d_forcing, sizeof(double) * (n_grids * n_species), cudaMemcpyDeviceToHost);

      // clean up
      cudaFree(d_rate_constants);
      cudaFree(d_state_variables);
      cudaFree(d_forcing);
      cudaFree(d_number_of_reactants);
      cudaFree(d_reactant_ids);
      cudaFree(d_number_of_products);
      cudaFree(d_product_ids);
      cudaFree(d_yields);
      cudaFree(device);
      return kernel_duration;
    }  // end of AddForcingTerms_kernelSetup
  }    // namespace cuda
}  // namespace micm