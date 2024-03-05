// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <chrono>
#include <iostream>
#include <micm/util/cuda_param.hpp>
#include <vector>

namespace micm
{
  namespace cuda
  {
    /// CUDA kernel to compute alpha - J[i] for each element i at the diagnoal of matrix J
    __global__ void AlphaMinusJacobianKernel(
        size_t n_grids,
        double* d_jacobian,
        size_t* d_jacobian_diagonal_elements,
        size_t jacobian_diagonal_elements_size,
        double alpha)
    {
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      if (tid < n_grids)
      {
        for (int j = 0; j < jacobian_diagonal_elements_size; j++)
        {
          size_t jacobian_index = d_jacobian_diagonal_elements[j];
          d_jacobian[jacobian_index + tid] += alpha;
        }
      }
    }

    // Specific CUDA device function to do reduction within a warp 
    // Use volatile to prevent compiler optimization (caching in registers)
    // No need to synchronize threads in the same warp
    __device__ void warpReduce(volatile double *sdata, size_t tid)
    {
      sdata[tid] += sdata[tid + 32];
      sdata[tid] += sdata[tid + 16];
      sdata[tid] += sdata[tid + 8];
      sdata[tid] += sdata[tid + 4];
      sdata[tid] += sdata[tid + 2];
      sdata[tid] += sdata[tid + 1];
    }

    // CUDA kernel to compute the scaled norm of the vector errors
    // Modified version from NVIDIA's reduction example: https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
    __global__ void NormalizedErrorKernel(
        double* d_y_old,
        double* d_y_new,
        double* d_errors,
        const size_t num_elements,
        const double atol,
        const double rtol,
        bool is_first_call)
    {
      // Declares a dynamically-sized shared memory array.
      // The size of this array is determined at runtime when the kernel is launched.
      // Shared memory is shared among all threads within the same block.
      extern __shared__ double sdata[];
    
      // Local thread ID within a threadblock
      size_t l_tid = threadIdx.x;
  
      // Global thread ID
      size_t g_tid = blockIdx.x*(blockDim.x*2) + threadIdx.x;

      if (is_first_call)
      {
        // Temporary device variables
        double d_ymax, d_scale;

        // Load two elements by one thread and do first add of reduction
        sdata[l_tid] = 0.0;
        for (int i=0; i<2; ++i)
        {
          if (g_tid < num_elements)
          {
            d_ymax = max(fabs(d_y_old[g_tid]), fabs(d_y_new[g_tid]));
            d_scale = atol + rtol * d_ymax;
            d_errors[g_tid] = d_errors[g_tid] * d_errors[g_tid] / (d_scale * d_scale);
            sdata[l_tid] += d_errors[g_tid];
          }
          g_tid += blockDim.x;
        } 
        __syncthreads();
      }
      else
      {
        // Load two elements by one thread and do first add of reduction
        // Access the d_errors array directly if it is not the first call
        sdata[l_tid] = 0.0;
        if (g_tid < num_elements) sdata[l_tid] += d_errors[g_tid];
        g_tid += blockDim.x;
        if (g_tid < num_elements) sdata[l_tid] += d_errors[g_tid];
        __syncthreads();
      }

      // Start at 1/2 block stride, do the add, and divide by two each iteration
      if (blockDim.x >= 2048) { if (l_tid < 1024) { sdata[l_tid] += sdata[l_tid + 1024]; } __syncthreads(); }
      if (blockDim.x >= 1024) { if (l_tid < 512)  { sdata[l_tid] += sdata[l_tid + 512];  } __syncthreads(); }
      if (blockDim.x >= 512)  { if (l_tid < 256)  { sdata[l_tid] += sdata[l_tid + 256];  } __syncthreads(); }
      if (blockDim.x >= 256)  { if (l_tid < 128)  { sdata[l_tid] += sdata[l_tid + 128];  } __syncthreads(); }
      if (blockDim.x >= 128)  { if (l_tid < 64)   { sdata[l_tid] += sdata[l_tid + 64];   } __syncthreads(); }
      if (l_tid < 32) warpReduce(sdata, l_tid);
  
      // Let the thread 0 of this threadblock write its result to output array, inexed by this threadblock
      if (l_tid == 0) d_errors[blockIdx.x] = sdata[0];
    }

    std::chrono::nanoseconds AlphaMinusJacobianDriver(
        CudaSparseMatrixParam& sparseMatrix,
        const std::vector<size_t> jacobian_diagonal_elements,
        double alpha)
    {
      // device pointers
      double* d_jacobian;
      size_t* d_jacobian_diagonal_elements;
      cudaMalloc(&d_jacobian, sizeof(double) * sparseMatrix.jacobian_size_);
      cudaMalloc(&d_jacobian_diagonal_elements, sizeof(size_t) * jacobian_diagonal_elements.size());
      cudaMemcpy(d_jacobian, sparseMatrix.jacobian_, sizeof(double) * sparseMatrix.jacobian_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_jacobian_diagonal_elements,
          jacobian_diagonal_elements.data(),
          sizeof(size_t) * jacobian_diagonal_elements.size(),
          cudaMemcpyHostToDevice);

      // kernel call
      size_t num_block = (sparseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      auto startTime = std::chrono::high_resolution_clock::now();
      AlphaMinusJacobianKernel<<<num_block, BLOCK_SIZE>>>(
          sparseMatrix.n_grids_, d_jacobian, d_jacobian_diagonal_elements, jacobian_diagonal_elements.size(), alpha);

      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);
      cudaMemcpy(sparseMatrix.jacobian_, d_jacobian, sizeof(double) * sparseMatrix.jacobian_size_, cudaMemcpyDeviceToHost);
      cudaFree(d_jacobian);
      cudaFree(d_jacobian_diagonal_elements);

      return kernel_duration;
    }

    double NormalizedErrorDriver(double* d_y_old, double* d_y_new, 
                                 double* d_errors, const size_t num_elements,
                                 const double atol, const double rtol)
    {
      size_t num_blocks = std::ceil(std::ceil(num_elements*1.0/BLOCK_SIZE) / 2.0);
      num_blocks = num_blocks < 1 ? 1 : num_blocks;
      size_t new_blocks;
      bool is_first_call;

      // kernel call
      // The "d_errors" vector will be overwritten by the reduction operation; need another temporary device vector if this is not desired
      is_first_call = true;
      NormalizedErrorKernel<<<num_blocks, BLOCK_SIZE, BLOCK_SIZE*sizeof(double)>>>(d_y_old, d_y_new, d_errors, 
                                                                                   num_elements, atol, rtol, is_first_call);

      is_first_call = false;
      while (num_blocks > 1)
      {
        // Update grid size
        new_blocks = std::ceil(std::ceil(num_blocks*1.0/BLOCK_SIZE) / 2.0);
        if (new_blocks <= 1)
        {
            NormalizedErrorKernel<<<1, BLOCK_SIZE, BLOCK_SIZE*sizeof(double)>>>(d_y_old, d_y_new, d_errors, 
                                                                                num_elements, atol, rtol,
                                                                                is_first_call);
            break;
        }
        NormalizedErrorKernel<<<new_blocks, BLOCK_SIZE, BLOCK_SIZE*sizeof(double)>>>(d_y_old, d_y_new, d_errors, 
                                                                                     num_elements, atol, rtol,
                                                                                     is_first_call);
        num_blocks = new_blocks;
      }
      cudaDeviceSynchronize();

      double error;
      cudaMemcpy(&error, &d_errors[0], sizeof(double), cudaMemcpyDeviceToHost);

      return std::max(std::sqrt(error / num_elements), 1.0e-10);
    }
  }  // namespace cuda
}  // namespace micm
