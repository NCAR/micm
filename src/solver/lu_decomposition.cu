// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#include <iostream> 
#include <vector>
#include <micm/util/cuda_param.hpp> 
struct DecomposeDevice{
    double* A_; 
    double* L_; 
    double* U_; 
    char* do_aik_; 
    size_t* aik_; 
    char* do_aki_;
    size_t* aki_;  
    size_t* uii_; 
    std::pair<size_t,size_t>* niLU_;
    std::pair<size_t, size_t>* uik_nkj_; 
    std::pair<size_t, size_t>* lij_ujk_;
    std::pair<size_t, size_t>* lki_nkj_; 
    std::pair<size_t, size_t>* lkj_uji_;
    size_t n_grids_;
    size_t niLU_size_;
}; 
namespace micm{
    namespace cuda{        
        __global__ void DecomposeKernel(DecomposeDevice* device)
        {
            size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
            double* A = device->A_; 
            double* L = device->L_;
            double* U = device->U_;
            std::pair<size_t, size_t>* lkj_uji = device->lkj_uji_;
            std::pair<size_t, size_t>* uik_nkj = device->uik_nkj_;
            std::pair<size_t, size_t>* lij_ujk = device->lij_ujk_;
            std::pair<size_t, size_t>* lki_nkj = device->lki_nkj_;
            size_t do_aik_offset = 0; //boolean vector 
            size_t aik_offset = 0;
            size_t uik_nkj_offset = 0; 
            size_t lij_ujk_offset = 0; 
            size_t do_aki_offset = 0; //boolean vector 
            size_t aki_offset = 0; 
            size_t lki_nkj_offset = 0; 
            size_t lkj_uji_offset = 0; 
            size_t uii_offset = 0; 
            
            if (tid < device->n_grids_){
                //loop through every element in niLU 
                for (size_t i = 0; i < device->niLU_size_; i++){
                    //upper triangular matrix 
                    auto inLU = device->niLU_[i]; 
                    for (size_t iU = 0; iU < inLU.second; ++iU){
                        if(device->do_aik_[do_aik_offset++]){
                            size_t U_idx = uik_nkj[uik_nkj_offset].first + tid;
                            size_t A_idx =  device->aik_[aik_offset++]+ tid; 
                            U[U_idx] = A[A_idx];
                        }
                        
                        for (size_t ikj = 0; ikj < uik_nkj[uik_nkj_offset].second; ++ikj){
                            size_t U_idx_1 = uik_nkj[uik_nkj_offset].first + tid; 
                            size_t L_idx = lij_ujk[lij_ujk_offset].first + tid;
                            size_t U_idx_2 = lij_ujk[lij_ujk_offset].second + tid;
                            U[U_idx_1] -= L[L_idx] * U[U_idx_2]; 
                           
                            ++lij_ujk_offset; 
                        }
                        ++uik_nkj_offset; 
                    }
                   // lower triangular matrix
                   
                    L[lki_nkj[lki_nkj_offset++].first + tid] = 1.0;         
                    
                    for (size_t iL = 0; iL <inLU.first; ++iL){
                        if(device->do_aki_[do_aki_offset++]){
                            size_t L_idx = lki_nkj[lki_nkj_offset].first + tid; 
                            size_t A_idx = device->aki_[aki_offset++] + tid; 
                            L[L_idx] = A[A_idx]; 
                        }
                        for(size_t ikj = 0; ikj < lki_nkj[lki_nkj_offset].second;++ikj){
                            size_t L_idx_1 = lki_nkj[lki_nkj_offset].first + tid;
                            size_t L_idx_2 = lkj_uji[lkj_uji_offset].first + tid;
                            size_t U_idx = lkj_uji[lkj_uji_offset].second + tid; 
                            L[L_idx_1] -= L[L_idx_2] * U[U_idx];
                            ++lkj_uji_offset; 
                        }
                        L[lki_nkj[lki_nkj_offset].first + tid]/=U[device->uii_[uii_offset] + tid]; 
                        ++lki_nkj_offset; 
                        ++uii_offset; 
                     }
                }
            }
        }// end of kernel
    
        void DecomposeKernelDriver(
            CudaSparseMatrixParam& sparseMatrix, 
            CudaSolverParam& solver){
            //create device pointers and allocate device memory 
            double* d_A; 
            double* d_L; 
            double* d_U; 
            bool* d_do_aik; 
            size_t* d_aik; 
            bool* d_do_aki;
            size_t* d_aki;  
            size_t* d_uii; 
            std::pair<size_t, size_t>* d_niLU; 
            std::pair<size_t, size_t>* d_uik_nkj; 
            std::pair<size_t, size_t>* d_lij_ujk;
            std::pair<size_t, size_t>* d_lki_nkj; 
            std::pair<size_t, size_t>* d_lkj_uji;
            DecomposeDevice* device; 
        
            cudaMalloc(&d_A,sizeof(double)* sparseMatrix.A_size_); 
            cudaMalloc(&d_L,sizeof(double)* sparseMatrix.L_size_); 
            cudaMalloc(&d_U,sizeof(double)* sparseMatrix.U_size_); 
            cudaMalloc(&d_do_aik,sizeof(char)* solver.do_aik_size_); 
            cudaMalloc(&d_aik,sizeof(size_t)* solver.aik_size_); 
            cudaMalloc(&d_do_aki,sizeof(char)* solver.do_aki_size_); 
            cudaMalloc(&d_aki,sizeof(size_t)* solver.aki_size_); 
            cudaMalloc(&d_uii,sizeof(size_t)* solver.uii_size_); 
            cudaMalloc(&d_niLU,sizeof(std::pair<size_t, size_t>)* solver.niLU_size_); 
            cudaMalloc(&d_uik_nkj,sizeof(std::pair<size_t, size_t>)* solver.uik_nkj_size_); 
            cudaMalloc(&d_lij_ujk,sizeof(std::pair<size_t, size_t>)* solver.lij_ujk_size_); 
            cudaMalloc(&d_lki_nkj,sizeof(std::pair<size_t, size_t>)* solver.lki_nkj_size_); 
            cudaMalloc(&d_lkj_uji,sizeof(std::pair<size_t, size_t>)* solver.lkj_uji_size_);
            cudaMalloc(&device, sizeof(DecomposeDevice)); 

            //transfer data from host to device 
            cudaMemcpy(d_A, sparseMatrix.A_, sizeof(double)* sparseMatrix.A_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_L, sparseMatrix.L_, sizeof(double)* sparseMatrix.L_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_U, sparseMatrix.U_, sizeof(double)* sparseMatrix.U_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aik, solver.do_aik_, sizeof(char)* solver.do_aik_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_aik, solver.aik_, sizeof(size_t)* solver.aik_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aki, solver.do_aki_, sizeof(char)* solver.do_aki_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_aki, solver.aki_, sizeof(size_t)*solver.aki_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_uii, solver.uii_, sizeof(size_t)* solver.uii_size_, cudaMemcpyHostToDevice);       
            cudaMemcpy(d_niLU, solver.niLU_, sizeof(std::pair<size_t, size_t>)*solver.niLU_size_, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_uik_nkj, solver.uik_nkj_, sizeof(std::pair<size_t, size_t>)*solver.uik_nkj_size_, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lij_ujk, solver.lij_ujk_, sizeof(std::pair<size_t, size_t>)*solver.lij_ujk_size_, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lki_nkj, solver.lki_nkj_, sizeof(std::pair<size_t, size_t>)*solver.lki_nkj_size_, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lkj_uji, solver.lkj_uji_, sizeof(std::pair<size_t, size_t>)*solver.lkj_uji_size_, cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->A_),&d_A, sizeof(double*), cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->L_),&d_L, sizeof(double*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->U_),&d_U, sizeof(double*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->do_aik_), &d_do_aik, sizeof(char*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->aik_), &d_aik, sizeof(size_t*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->do_aki_),&d_do_aki,sizeof(char*),cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->aki_),&d_aki, sizeof(size_t*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->uii_), &d_uii, sizeof(size_t*), cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->niLU_), &d_niLU, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->uik_nkj_), &d_uik_nkj, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->lij_ujk_), &d_lij_ujk, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->lki_nkj_), &d_lki_nkj, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->lkj_uji_), &d_lkj_uji, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            
            //total number of threads is number of blocks in sparseMatrix A 
            size_t num_block = (sparseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE; 
            device->n_grids_ = sparseMatrix.n_grids_;  
            device->niLU_size_ = solver.niLU_size_; 

           // call kernel
            DecomposeKernel<<<num_block, BLOCK_SIZE>>>(device); 
            cudaDeviceSynchronize();
            cudaMemcpy(sparseMatrix.L_, d_L, sizeof(double)* sparseMatrix.L_size_, cudaMemcpyDeviceToHost); 
            cudaMemcpy(sparseMatrix.U_, d_U, sizeof(double)* sparseMatrix.U_size_, cudaMemcpyDeviceToHost); 
          
        //clean up 
        cudaFree(d_A); 
        cudaFree(d_L); 
        cudaFree(d_U); 
        cudaFree(d_do_aik); 
        cudaFree(d_aik);
        cudaFree(d_do_aki); 
        cudaFree(d_aki); 
        cudaFree(d_uii); 
        cudaFree(device); 
    }//end kernelDriver
 }//end cuda 
}//end micm