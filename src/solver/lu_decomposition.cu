#include <iostream> 
#include <micm/util/cuda_param> 
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
        __global__ void Decompose_kernel(){

        }
    
        void DecomposeKernelDriver(
            CUDAMatrixParam& matrix, 
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

            cudaMalloc(&d_A, sizeof(double)* matrix.A_size); 
            cudaMalloc(&d_L, sizeof(double)* matrix.L_size); 
            cudaMalloc(&d_U, sizeof(double)* matrix.U_size); 
            cudaMalloc(&d_do_aik, sizeof(bool)* solver.do_aik_size); 
            cudaMalloc(&d_aik, sizeof(size_t)* solver.aik_size); 
            cudaMalloc(&d_do_aki, sizeof(bool)* solver.do_aki_size); 
            cudaMalloc(&d_aki, sizeof(size_t)* solver.aki_size); 
            cudaMalloc(&d_uii, sizeof(size_t)* solver.uii_size); 
            cudaMalloc(&device, sizeof(decomposeDevice)); 

            //transfer data from host to device 
            cudaMemcpy(d_A, matrix.A, sizeof(double)* matrix.A_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_L, matrix.L, sizeof(double)* matrix.L_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_U, matrix.U, sizeof(double)* matrix.U_size, cudaMemcpyHostToDevice); 
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
            
            
        

            }
        


    }
}