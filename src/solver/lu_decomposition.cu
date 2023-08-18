#include <iostream> 
struct pair{
    size_t first; 
    size_t second; 
}
namespace micm{
    namespace cuda{
        __global__ void Decompose_kernel(){

        }
    
        void decompose_kernelSetup(
            const double* A, 
            size_t A_size,
            double* L, 
            size_t L_size, 
            double* U,
            size_t U_size,
            size_t* niLu, //pair 
            size_t* niLu_size,
            bool* do_aik,
            size_t do_aik_size
            size_t* aik, 
            size_t aik_size,
            size_t* uik_nkj, //pair 
            size_t uik_nki_size,
            size_t* lij_ujk, //pair 
            size_t lij_ujk_size,
            bool* do_aki, 
            size_t do_aki_size,
            size_t* aki, 
            size_t aki_size,
            size_t* lki_nkj, //pair 
            size_t lki_nkj_size,
            size_t* lkj_uji, //pair
            size_t lki_uji_size,
            size_t* uii,
            size_t uii_size){
            
            //create device pointers and allocate device memory 
            double* d_A; 
            double* d_L; 
            double* d_U; 
            bool* d_do_aik; 
            size_t* d_aik; 
            bool* d_do_aki;
            size_t* d_aki;  
            size_t* d_uii; 

            cudaMalloc(&d_A, A, sizeof(double)* A_size); 
            cudaMalloc(&d_L, L, sizeof(double)* L_size); 
            cudaMalloc(&d_U, U, sizeof(double)* U_size); 
            cudaMalloc(&d_do_aik, sizeof(bool)* do_aik_size); 
            cudaMalloc(&d_aik, sizeof(size_t)* aik_size); 
            cudaMalloc(&d_do_aki, sizeof(bool)* do_aki_size); 
            cudaMalloc(&d_aki, sizeof(size_t)* aki_size); 
            cudaMalloc(d_uii, sizeof(size_t)* uii_size); 

            //transfer data from host to device 
            cudaMemcpy(d_A, A, sizeof(double)* A_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_L, L, sizeof(double)* L_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_U, U, sizeof(double)* U_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aik, do_aik, sizeof(bool)* do_aik_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_aik, aik, sizeof(size_t)* aik_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aki, do_aki, sizeof(bool)* do_aki_size, cudaMemcpyHostToDevice); 

            


            }
        


    }
}