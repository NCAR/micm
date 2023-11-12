// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#include <iostream> 
#include <vector>
#include <chrono>
#include <micm/util/cuda_param.hpp> 
struct SolveDevice{
    std::pair<size_t, size_t>* nLij_Lii_;
    std::pair<size_t, size_t>* Lij_yj_; 
    std::pair<size_t, size_t>* nUij_Uii_;
    std::pair<size_t, size_t>* Uij_xj_;
    double* lower_matrix_;
    double* upper_matrix_; 
    double* b_; 
    double* x_;
};
namespace micm{
    namespace cuda{
__global__ void SolveKernel(SolveDevice* device, 
                            size_t n_grids,
                            size_t b_column_counts,
                            size_t x_column_counts,
                            size_t nLij_Lii_size,
                            size_t nUij_Uii_size)
{
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    double* b = device->b_;
    double* x = device->x_;
    double* y = device->x_; //Alias x for consistency with equation, but to reuse memory
    double* lower_matrix = device->lower_matrix_;
    double* upper_matrix = device->upper_matrix_;
    std::pair<size_t, size_t>* nLij_Lii = device->nLij_Lii_;
    std::pair<size_t, size_t>* Lij_yj = device->Lij_yj_;
    std::pair<size_t, size_t>* nUij_Uii = device->nUij_Uii_;
    std::pair<size_t, size_t>* Uij_xj = device->Uij_xj_;

   if (tid < n_grids)
   { 
        size_t b_column_index = 0;
        size_t x_column_index = 0;
        size_t y_column_index = 0;
        size_t b_column_backward_index = b_column_counts -1;
        size_t x_column_backward_index = x_column_counts -1;
        size_t Lij_yj_index = 0; 
        size_t Uij_xj_index = 0;
       
        for (size_t j = 0; j < nLij_Lii_size; ++j)
        {
            auto& nLij_Lii_element = nLij_Lii[j]; 
            y[y_column_index * n_grids + tid] = b[b_column_index++ * n_grids + tid]; 
            for (size_t i = 0; i < nLij_Lii_element.first; ++i)
            {
                size_t lower_matrix_index = Lij_yj[Lij_yj_index].first + tid;
                size_t y_index = Lij_yj[Lij_yj_index].second * n_grids + tid; 
                y[y_column_index * n_grids + tid] -= lower_matrix[lower_matrix_index] * y[y_index];
                ++Lij_yj_index;  
            }
            y[y_column_index++ * n_grids + tid] /= lower_matrix[nLij_Lii_element.second + tid]; 
        }
        
        for (size_t k = 0; k < nUij_Uii_size; ++k)
        {   
            auto& nUij_Uii_element = nUij_Uii[k]; 
        
            for (size_t i = 0; i < nUij_Uii_element.first; ++i)
            {
                size_t upper_matrix_index = Uij_xj[Uij_xj_index].first + tid;
                size_t x_index = Uij_xj[Uij_xj_index].second * n_grids + tid;
                x[x_column_backward_index * n_grids + tid] -= upper_matrix[upper_matrix_index] * x[x_index];
                ++Uij_xj_index;
            }
            x[x_column_backward_index * n_grids + tid] /= upper_matrix[nUij_Uii_element.second + tid];
            
            if (x_column_backward_index != 0)
            {
                --x_column_backward_index;
            }
        }
    }
}
    std::chrono::nanoseconds SolveKernelDriver(CudaLinearSolverParam& linearSolver,CudaSparseMatrixParam& sparseMatrix, CudaMatrixParam& denseMatrix)
    {
    //create device pointer
    std::pair<size_t, size_t>* d_nLij_Lii; 
    std::pair<size_t, size_t>* d_Lij_yj; 
    std::pair<size_t, size_t>* d_nUij_Uii; 
    std::pair<size_t, size_t>* d_Uij_xj;
    double* d_lower_matrix; 
    double* d_upper_matrix;
    double* d_b; 
    double* d_x;
    SolveDevice* device;

    //allocate device memory 
    cudaMalloc(&d_nLij_Lii, sizeof(std::pair<size_t, size_t>)* linearSolver.nLij_Lii_size_); 
    cudaMalloc(&d_Lij_yj, sizeof(std::pair<size_t, size_t>)* linearSolver.Lij_yj_size_); 
    cudaMalloc(&d_nUij_Uii, sizeof(std::pair<size_t, size_t>)* linearSolver.nUij_Uii_size_);
    cudaMalloc(&d_Uij_xj, sizeof(std::pair<size_t, size_t>)* linearSolver.Uij_xj_size_); 
    
    cudaMalloc(&d_lower_matrix, sizeof(double)* sparseMatrix.lower_matrix_size_); 
    cudaMalloc(&d_upper_matrix, sizeof(double)* sparseMatrix.upper_matrix_size_);
    cudaMalloc(&d_b, sizeof(double)* denseMatrix.b_size_);
    cudaMalloc(&d_x, sizeof(double)* denseMatrix.x_size_); 
    cudaMalloc(&device, sizeof(SolveDevice));

    //transfer memory from host to device
    cudaMemcpy(d_nLij_Lii, linearSolver.nLij_Lii_, sizeof(std::pair<size_t, size_t>)* linearSolver.nLij_Lii_size_,cudaMemcpyHostToDevice);
    cudaMemcpy(d_Lij_yj, linearSolver.Lij_yj_, sizeof(std::pair<size_t, size_t>)* linearSolver.Lij_yj_size_,cudaMemcpyHostToDevice);
    cudaMemcpy(d_nUij_Uii, linearSolver.nUij_Uii_, sizeof(std::pair<size_t, size_t>)* linearSolver.nUij_Uii_size_,cudaMemcpyHostToDevice);
    cudaMemcpy(d_Uij_xj, linearSolver.Uij_xj_, sizeof(std::pair<size_t, size_t>)* linearSolver.Uij_xj_size_, cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_lower_matrix, sparseMatrix.lower_matrix_, sizeof(double)*sparseMatrix.lower_matrix_size_, cudaMemcpyHostToDevice);
    cudaMemcpy(d_upper_matrix, sparseMatrix.upper_matrix_, sizeof(double)*sparseMatrix.upper_matrix_size_, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, denseMatrix.b_, sizeof(double)* denseMatrix.b_size_, cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, denseMatrix.x_, sizeof(double)* denseMatrix.x_size_, cudaMemcpyHostToDevice);
    
    cudaMemcpy(&(device->nLij_Lii_), &d_nLij_Lii, sizeof(std::pair<size_t, size_t>*),cudaMemcpyHostToDevice);
    cudaMemcpy(&(device->Lij_yj_), &d_Lij_yj, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(device->nUij_Uii_), &d_nUij_Uii, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(device->Uij_xj_), &d_Uij_xj, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice);
    
    cudaMemcpy(&(device->lower_matrix_), &d_lower_matrix, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(device->upper_matrix_), &d_upper_matrix, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(&(device->b_), &d_b, sizeof(double*), cudaMemcpyHostToDevice); 
    cudaMemcpy(&(device->x_),&d_x, sizeof(double*), cudaMemcpyHostToDevice);
    
    //kernel call 
    size_t num_block = (denseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
    auto startTime = std::chrono::high_resolution_clock::now();
    SolveKernel<<<num_block, BLOCK_SIZE>>>(device, 
                                          denseMatrix.n_grids_, 
                                          denseMatrix.b_column_counts_, 
                                          denseMatrix.x_column_counts_,
                                          linearSolver.nLij_Lii_size_,
                                          linearSolver.nUij_Uii_size_);
    cudaDeviceSynchronize();
    auto endTime = std::chrono::high_resolution_clock::now();
    auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);
    cudaMemcpy(denseMatrix.x_, d_x, sizeof(double)* denseMatrix.x_size_, cudaMemcpyDeviceToHost);

    //clean up 
    cudaFree(d_nLij_Lii);
    cudaFree(d_Lij_yj); 
    cudaFree(d_nUij_Uii); 
    cudaFree(d_Uij_xj);
    cudaFree(d_lower_matrix); 
    cudaFree(d_upper_matrix);
    cudaFree(d_b); 
    cudaFree(d_x);
    cudaFree(device);
    return kernel_duration;
    }
  }//end cuda 
}// end micm 