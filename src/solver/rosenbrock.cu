#pragma once
#include <chrono>
#include <vector>
#include <iostream>
#include <micm/util/cuda_param.hpp>

namespace micm{
    namespace cuda{
        __global__ void AlphaMinusJacobianKernel(size_t n_grids,
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
        
        std::chrono::nanoseconds AlphaMinusJacobianDriver(
                        CudaSparseMatrixParam& sparseMatrix,
                        const std::vector<size_t> jacobian_diagonal_elements,
                        double alpha)
    {
        //device pointers
        double* d_jacobian;
        size_t* d_jacobian_diagonal_elements; 
        cudaMalloc(&d_jacobian, sizeof(double)* sparseMatrix.jacobian_size_); 
        cudaMalloc(&d_jacobian_diagonal_elements, sizeof(size_t)* jacobian_diagonal_elements.size());
        cudaMemcpy(d_jacobian, sparseMatrix.jacobian_, sizeof(double)* sparseMatrix.jacobian_size_, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_jacobian_diagonal_elements, jacobian_diagonal_elements.data(), sizeof(size_t)* jacobian_diagonal_elements.size(), cudaMemcpyHostToDevice);

        //kernel call
        size_t num_block = (sparseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        auto startTime = std::chrono::high_resolution_clock::now();
        AlphaMinusJacobianKernel<<<num_block, BLOCK_SIZE>>>(sparseMatrix.n_grids_,
                                d_jacobian,  
                                d_jacobian_diagonal_elements,
                                jacobian_diagonal_elements.size(),
                                alpha);


        cudaDeviceSynchronize();
        auto endTime = std::chrono::high_resolution_clock::now();
        auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);
        cudaMemcpy(sparseMatrix.jacobian_, d_jacobian, sizeof(double)* sparseMatrix.jacobian_size_, cudaMemcpyDeviceToHost);
        cudaFree(d_jacobian);
        cudaFree(d_jacobian_diagonal_elements);
        
        return kernel_duration; 

    }

    }// end cuda
}// end mimc 
