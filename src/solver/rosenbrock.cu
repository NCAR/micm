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
            //printf("j: %d\n", j); 
            size_t jacobian_index = d_jacobian_diagonal_elements[j];
            printf("jacobian index: %d\n", d_jacobian_diagonal_elements[j]); 
            d_jacobian[jacobian_index + tid] += alpha; 
        }
    } 
}
        
        void AlphaMinusJacobianDriver(
                        CudaSparseMatrixParam& sparseMatrix,
                        
                        double alpha)
    {
        //device pointers
        std::cout<< "element size: "<<sparseMatrix.jacobian_diagonal_elements_size_<<std::endl; 
        
        
        double* d_jacobian;
        size_t* d_jacobian_diagonal_elements; 
        cudaMalloc(&d_jacobian, sizeof(double)* sparseMatrix.jacobian_size_); 
        cudaMalloc(&d_jacobian_diagonal_elements, sizeof(size_t)*sparseMatrix.jacobian_diagonal_elements_size_);
        cudaMemcpy(d_jacobian, sparseMatrix.jacobian_, sparseMatrix.jacobian_size_, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_jacobian_diagonal_elements, sparseMatrix.jacobian_diagonal_elements_, sparseMatrix.jacobian_diagonal_elements_size_, cudaMemcpyHostToDevice);
        
        
        //kernel call
        size_t num_block = (sparseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
        AlphaMinusJacobianKernel<<<num_block, BLOCK_SIZE>>>(sparseMatrix.n_grids_,
                                d_jacobian,  
                                d_jacobian_diagonal_elements,
                                sparseMatrix.jacobian_diagonal_elements_size_,
                                alpha);
        
        cudaDeviceSynchronize();
        cudaMemcpy(sparseMatrix.jacobian_, d_jacobian, sparseMatrix.jacobian_size_, cudaMemcpyDeviceToHost);
        cudaFree(d_jacobian);
        cudaFree(d_jacobian_diagonal_elements);
    }

    }// end cuda
}// end mimc 
