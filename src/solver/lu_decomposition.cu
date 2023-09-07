#include <iostream> 
#include <micm/util/cuda_param> 
#include<thrust/device_vector.h> 
#include <thrust/pair.h>
const BLOCK_SIZE = 320; 
struct decomposeDevice{
    double* A; 
    double* L; 
    double* U; 
    bool* do_aik; 
    size_t* aik; 
    bool* do_aki;
    size_t* aki;  
    size_t* uii; 
}; 

namespace micm{
    namespace cuda{
        __global__ void DecomposeKernel(
            decomposeDevice& device, 
            thrust::device_vector<thrust::pair<size_t,size_t>> niLU,
            thrust::device_vector<thrust::pair<size_t,size_t>> uik_nkj, 
            thrust::device_vector<thrust::pair<size_t,size_t>> lij_ujk,
            thrust::device_vector<thrust::pair<size_t,size_t>> lki_nkj,
            thrust::device_vector<thrust::pair<size_t,size_t>> lkj_uji)
        {
            size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
            double* A = device->A; 
            double* L = device->L;
            double* U = device->U;
            size_t do_aik_offset = 0; //boolean vector 
            size_t aik_offset = 0;
            size_t uik_nkj_offset = 0; 
            size_t lij_ujk_offset = 0; 
            size_t do_aki_offset = 0; //boolean vector 
            size_t aki_offset = 0; 
            size_t lki_nkj_offset = 0; 
            size_t lkj_uji_offset = 0; 
            size_t uii_offset = 0; 
            if (tid < A_size){
                for (auto& inLU : niLU){
                    //upper triangular matrix 
                    for (size_t iU = 0; iU < inLU.second; ++iU){
                        if(device->do_aik[++do_aik_offset]){
                            size_t U_idx = uik_nkj[uik_nkj_offset]->first + tid;
                            size_t A_idx =  device->aik[++aik_offset]+ tid; 
                            U[U_idx] = A[A_idx]; 
                        }
                        for (size_t ikj = 0; ikj < uik_nkj[uik_nkj_offset]->second; ++ikj){
                            
                            size_t L_idx = lij_ujk[lij_ujk_offset]->first + tid;
                            size_t U_idx_1 = uik_nkj[uik_nkj_offset]->first + tid; 
                            size_t U_idx_2 = lij_ujk[lij_ujk_offset]->second + tid; 
                            U[U_idx_1] -= L[L_idx] * U[U_idx_2]; 
                            ++lij_ujk_offset; 
                        }
                        ++uik_nkj_offset; 
                    }
                    //lower triangular matrix
                    L[lki_nkj[++lki_nkj_offset]->first + tid] = 1.0; 
                    for (size_t iL = 0; iL <inLU.first; ++iL){
                        if(device->do_aki[++do_aki_offset]){
                            size_t L_idx = lki_nkj[lkj_nkj_offset]->first + tid; 
                            size_t A_idx = aki->device[++aki_offset] + tid; 
                            L[L_idx] = A[A_idx]; 
                        }
                        //working in progress 
                        for(size_t ikj = 0; ikj < lki_nkj[lki_nkj_offset]->second;++ikj){
                            size_t L_idx_1 = lki_nkj[lki_nkj_offset]->first + tid;
                            size_t L_idx_2 = lkj_uji[lkj_uji_offset]->first + tid;
                            size_t U_idx = lkj_uji[lkj_uji_offset]->second + tid; 
                            ++lkj_uji_offset; 
                        }
                        size_t L_idx = lki_nkj[lki_nkj_offset]->first + tid; 
                        size_t U_idx = device->uii[uii_offset]+tid; 
                        L[L_idx]/=U[U_idx]; 
                        ++lki_nkj_offset; 
                        ++uii_offset; 
                    }
                }
            }
        }// end of kernel
    
        void DecomposeKernelDriver(
            CUDAMatrixParam& sparseMatrix, 
            CUDASolverParam& solver){
            
            //create device pointers and allocate device memory 
            double* d_A; 
            double* d_L; 
            double* d_U; 
            bool* d_do_aik; 
            size_t* d_aik; 
            bool* d_do_aki;
            size_t* d_aki;  
            size_t* d_uii; 
            decomposeDevice* device; 

            cudaMalloc(&d_A, sizeof(double)* sparseMatrix.A_size); 
            cudaMalloc(&d_L, sizeof(double)* sparseMatrix.L_size); 
            cudaMalloc(&d_U, sizeof(double)* sparseMatrix.U_size); 
            cudaMalloc(&d_do_aik, sizeof(bool)* solver.do_aik_size); 
            cudaMalloc(&d_aik, sizeof(size_t)* solver.aik_size); 
            cudaMalloc(&d_do_aki, sizeof(bool)* solver.do_aki_size); 
            cudaMalloc(&d_aki, sizeof(size_t)* solver.aki_size); 
            cudaMalloc(&d_uii, sizeof(size_t)* solver.uii_size); 
            cudaMalloc(&device, sizeof(decomposeDevice)); 

            //transfer data from host to device 
            cudaMemcpy(d_A, sparseMatrix.A, sizeof(double)* sparseMatrix.A_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_L, sparseMatrix.L, sizeof(double)* sparseMatrix.L_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_U, sparseMatrix.U, sizeof(double)* sparseMatrix.U_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aik, solver.do_aik, sizeof(bool)* solver.do_aik_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_aik, solver.aik, sizeof(size_t)* solver.aik_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aki, solver.do_aki, sizeof(bool)* solver.do_aki_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_uii, solver.uii, sizeof(size_t)* solver.uii_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->A),&d_A, sizeof(double*), cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->L),&d_L, sizeof(double*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->U),&d_U, sizeof(double*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->do_aik), &d_do_aik, sizeof(bool*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->aik), &d_aik, sizeof(size_t*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->do_aki),&d_do_aki,sizeof(bool*),cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->aki),&d_aki, sizeof(size_t*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->uii), &d_uii, sizeof(size_t*), cudaMemcpyHostToDevice);
            
            
            size_t num_block = (sparseMatrix.A_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
            //call kernel
            DecomposeKernel<<<BLOCK_SIZE, num_block>>>(
            decomposeDevice& device, 
            niLU,
            uik_nkj, 
            lij_ujk,
            lki_nkj,
            lkj_uji)

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