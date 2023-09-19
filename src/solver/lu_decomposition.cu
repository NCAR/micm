#include <iostream> 
#include <vector>
#include <micm/util/cuda_param.hpp> 
const size_t BLOCK_SIZE = 320; 
struct decomposeDevice{
    double* A; 
    double* L; 
    double* U; 
    char* do_aik; 
    size_t* aik; 
    char* do_aki;
    size_t* aki;  
    size_t* uii; 
    std::pair<size_t,size_t>* niLU;
    std::pair<size_t, size_t>* uik_nkj; 
    std::pair<size_t, size_t>* lij_ujk;
    std::pair<size_t, size_t>* lki_nkj; 
    std::pair<size_t, size_t>* lkj_uji;
}; 
namespace micm{
    namespace cuda{
        __global__ void pairCheck(decomposeDevice* device, size_t niLU_size){
            size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
            if (tid < niLU_size){
            auto niLU_element = device->niLU[tid].second;
            printf("The second element of niLU is %d\n", niLU_element);}//print something
        }
        
        
        __global__ void DecomposeKernel(
            decomposeDevice* device,
            size_t n_grids,
            size_t niLU_size)
        {
            size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
            double* A = device->A; 
            double* L = device->L;
            double* U = device->U;
            std::pair<size_t, size_t>* lkj_uji = device->lkj_uji;
            std::pair<size_t, size_t>* uik_nkj = device->uik_nkj;
            std::pair<size_t, size_t>* lij_ujk = device->lij_ujk;
            std::pair<size_t, size_t>* lki_nkj = device->lki_nkj;
            size_t do_aik_offset = 0; //boolean vector 
            size_t aik_offset = 0;
            size_t uik_nkj_offset = 0; 
            size_t lij_ujk_offset = 0; 
            size_t do_aki_offset = 0; //boolean vector 
            size_t aki_offset = 0; 
            size_t lki_nkj_offset = 0; 
            size_t lkj_uji_offset = 0; 
            size_t uii_offset = 0; 
            if (tid < n_grids){
        
                //loop through every element in niLU 
                for (size_t i = 0; i < niLU_size; i++){
                    //upper triangular matrix 
                    auto inLU = device->niLU[i]; 
                    for (size_t iU = 0; iU < inLU.second; ++iU){
                        if(device->do_aik[do_aik_offset++]){
                            prinf("tid: %d, this is aik\n",device->aik[aik_offset]);
                            size_t U_idx = uik_nkj[uik_nkj_offset].first + tid;
                            size_t A_idx =  device->aik[aik_offset++]+ tid; 
                            prinf("tid: %d, this is uik_nkj_first\n",uik_nkj[uik_nkj_offset].first);
                            printf("tid: %d, this is gpu u index: %d\n", tid, U_idx); 
                            printf("tid: %d, this is gpu A index: %d\n", tid, A_idx); 
                            U[U_idx] = A[A_idx]; 
                            printf ("tid: %d, this is gpu U value: %d\n",tid, U[U_idx]); 
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
                    //lower triangular matrix
                    L[lki_nkj[lki_nkj_offset++].first + tid] = 1.0; 
                    
                    for (size_t iL = 0; iL <inLU.first; ++iL){
                        if(device->do_aki[do_aki_offset++]){
                            size_t L_idx = lki_nkj[lki_nkj_offset].first + tid; 
                            size_t A_idx = device->aki[aki_offset++] + tid; 
                            L[L_idx] = A[A_idx]; 
                        }
                        //working in progress 
                        for(size_t ikj = 0; ikj < lki_nkj[lki_nkj_offset].second;++ikj){
                            size_t L_idx_1 = lki_nkj[lki_nkj_offset].first + tid;
                            size_t L_idx_2 = lkj_uji[lkj_uji_offset].first + tid;
                            size_t U_idx = lkj_uji[lkj_uji_offset].second + tid; 
                            L[L_idx_1] -= L[L_idx_2] * U[U_idx];
                            ++lkj_uji_offset; 
                        }
                        size_t L_idx = lki_nkj[lki_nkj_offset].first + tid; 
                        size_t U_idx = device->uii[uii_offset] + tid; 
                        L[L_idx]/=U[U_idx]; 
                        ++lki_nkj_offset; 
                        ++uii_offset; 
                    }
                }
            }
        }// end of kernel
    
        void DecomposeKernelDriver(
            CUDASparseMatrixParam& sparseMatrix, 
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
            std::pair<size_t, size_t>* d_niLU; 
            std::pair<size_t, size_t>* d_uik_nkj; 
            std::pair<size_t, size_t>* d_lij_ujk;
            std::pair<size_t, size_t>* d_lki_nkj; 
            std::pair<size_t, size_t>* d_lkj_uji;
            decomposeDevice* device; 
        
            cudaMalloc(&d_A,sizeof(double)* sparseMatrix.A_size); 
            cudaMalloc(&d_L,sizeof(double)* sparseMatrix.L_size); 
            cudaMalloc(&d_U,sizeof(double)* sparseMatrix.U_size); 
            cudaMalloc(&d_do_aik,sizeof(char)* solver.do_aik_size); 
            cudaMalloc(&d_aik,sizeof(size_t)* solver.aik_size); 
            cudaMalloc(&d_do_aki,sizeof(char)* solver.do_aki_size); 
            cudaMalloc(&d_aki,sizeof(size_t)* solver.aki_size); 
            cudaMalloc(&d_uii,sizeof(size_t)* solver.uii_size); //check upon here
            cudaMalloc(&d_niLU,sizeof(std::pair<size_t, size_t>)* solver.niLU_size); 
            cudaMalloc(&d_uik_nkj,sizeof(std::pair<size_t, size_t>)* solver.uik_nkj_size); 
            cudaMalloc(&d_lij_ujk,sizeof(std::pair<size_t, size_t>)* solver.lij_ujk_size); 
            cudaMalloc(&d_lki_nkj,sizeof(std::pair<size_t, size_t>)* solver.lki_nkj_size); 
            cudaMalloc(&d_lkj_uji,sizeof(std::pair<size_t, size_t>)* solver.lkj_uji_size);
            cudaMalloc(&device, sizeof(decomposeDevice)); 

            //transfer data from host to device 
            cudaMemcpy(d_A, sparseMatrix.A, sizeof(double)* sparseMatrix.A_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_L, sparseMatrix.L, sizeof(double)* sparseMatrix.L_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_U, sparseMatrix.U, sizeof(double)* sparseMatrix.U_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aik, solver.do_aik, sizeof(char)* solver.do_aik_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_aik, solver.aik, sizeof(size_t)* solver.aik_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_do_aki, solver.do_aki, sizeof(char)* solver.do_aki_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_uii, solver.uii, sizeof(size_t)* solver.uii_size, cudaMemcpyHostToDevice);       
            cudaMemcpy(d_niLU, solver.niLU, sizeof(std::pair<size_t, size_t>)*solver.niLU_size, cudaMemcpyHostToDevice); 
            cudaMemcpy(d_uik_nkj, solver.uik_nkj, sizeof(std::pair<size_t, size_t>)*solver.uik_nkj_size, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lij_ujk, solver.lij_ujk, sizeof(std::pair<size_t, size_t>)*solver.lij_ujk_size, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lki_nkj, solver.lki_nkj, sizeof(std::pair<size_t, size_t>)*solver.lki_nkj_size, cudaMemcpyHostToDevice);
            cudaMemcpy(d_lkj_uji, solver.lkj_uji, sizeof(std::pair<size_t, size_t>)*solver.lkj_uji_size, cudaMemcpyHostToDevice);

            cudaMemcpy(&(device->A),&d_A, sizeof(double*), cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->L),&d_L, sizeof(double*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->U),&d_U, sizeof(double*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->do_aik), &d_do_aik, sizeof(char*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->aik), &d_aik, sizeof(size_t*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->do_aki),&d_do_aki,sizeof(char*),cudaMemcpyHostToDevice); //check upon here 
            cudaMemcpy(&(device->aki),&d_aki, sizeof(size_t*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->uii), &d_uii, sizeof(size_t*), cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->niLU), &d_niLU, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice);
            cudaMemcpy(&(device->uik_nkj), &d_uik_nkj, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->lij_ujk), &d_lij_ujk, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->lki_nkj), &d_lki_nkj, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            cudaMemcpy(&(device->lkj_uji), &d_lkj_uji, sizeof(std::pair<size_t, size_t>*), cudaMemcpyHostToDevice); 
            
            //total number of threads is number of blocks in sparseMatrix A 
            size_t num_block = (sparseMatrix.n_grids + BLOCK_SIZE - 1) / BLOCK_SIZE; 
            size_t n_grids = sparseMatrix.n_grids;  
            size_t niLU_size = solver.niLU_size; 
         
            //call kernel
            DecomposeKernel<<<num_block, BLOCK_SIZE>>>(device, n_grids, niLU_size); 
            cudaDeviceSynchronize();
            cudaMemcpy(sparseMatrix.L, d_L, sizeof(double)* sparseMatrix.L_size, cudaMemcpyDeviceToHost); 
            cudaMemcpy(sparseMatrix.U, d_U, sizeof(double)* sparseMatrix.U_size, cudaMemcpyDeviceToHost); 
           
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