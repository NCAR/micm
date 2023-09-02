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
        // __global__ void Decompose_kernel(
        //     decomposeDevice& device, 
        //     thrust::device_vector d_niLU<thrust::pair<size_t,size_t>>;
        // )
        // {

        // }
    
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

            }
        


    }
}