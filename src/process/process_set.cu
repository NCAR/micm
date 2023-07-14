#include <iostream>

namespace micm {
    namespace cuda {
      // flipped memory layout
    __global__ void AddForcingTerms_kernel(
        double* rate_constants, 
        double* state_variables, 
        double* forcing, 
        int ngrids,
        int nrxns, 
        int nspecs, 
        size_t* number_of_reactants_, 
        size_t* reactant_ids_, 
        size_t* number_of_products_, 
        size_t* product_ids_, 
        double* yields_)
    {
      //define thread index 
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x; 
      size_t react_id_offset, prod_id_offset, yield_offset;

      if (tid < ngrids){
         react_id_offset = 0;
	 prod_id_offset = 0;
	 yield_offset = 0;
         for (std::size_t i_rxn = 0; i_rxn < nrxns; ++i_rxn)
         {
           double rate = rate_constants[i_rxn*ngrids+tid];
           size_t index;

           for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
             rate *= state_variables[reactant_ids_[react_id_offset+i_react]*ngrids+tid];
           for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react){
             forcing[reactant_ids_[react_id_offset+i_react]*ngrids+tid] -= rate;
	   }
           for (std::size_t i_prod = 0; i_prod < number_of_products_[i_rxn]; ++i_prod){
             index = product_ids_[prod_id_offset+i_prod]*ngrids+tid;
             forcing[index] += yields_[yield_offset+i_prod] * rate;
           }
           react_id_offset += number_of_reactants_[i_rxn];
           prod_id_offset += number_of_products_[i_rxn];
           yield_offset += number_of_products_[i_rxn];
         } // for loop over number of reactions
      }    // if check for valid CUDA threads
    }      // end of AddForcingTerms_kernel

    void AddForcingTerms_kernelSetup(
        const double* rate_constants_data,
        const double* state_variables_data,
        double* forcing_data,
        int ngrids,
        int nrxns,
        int nspecs,
        const size_t* number_of_reactants,
        int number_of_reactants_size,
        const size_t* reactant_ids, 
        int reactant_ids_size,
        const size_t* number_of_products, 
        int number_of_products_size,
        const size_t* product_ids,
        int product_ids_size,
        const double* yields,
        int yields_size){
        
        // device pointer to vectorss
        double* d_rate_constants; 
        double* d_state_variables; 
        double* d_forcing; 
        double* d_yields_; 
        size_t* d_number_of_reactants_; 
        size_t* d_reactant_ids_; 
        size_t* d_number_of_products_; 
        size_t* d_product_ids_; 
       
        //allocate device memory
        size_t rate_constants_bytes = sizeof(double) * (ngrids * nrxns); 
        size_t state_forcing_bytes = sizeof(double) * (ngrids * nspecs); 
        size_t yields_bytes = sizeof(double) * yields_size;
        size_t number_of_reactants_bytes = sizeof(size_t) * number_of_reactants_size;
        size_t reactant_ids_bytes = sizeof(size_t) * reactant_ids_size; 
        size_t number_of_products_bytes = sizeof(size_t) * number_of_products_size; 
        size_t product_ids_bytes = sizeof(size_t) * product_ids_size;
        
        cudaMalloc(&d_rate_constants, rate_constants_bytes); 
        cudaMalloc(&d_state_variables, state_forcing_bytes); 
        cudaMalloc(&d_forcing, state_forcing_bytes); 
        cudaMalloc(&d_number_of_reactants_, number_of_reactants_bytes);
        cudaMalloc(&d_reactant_ids_, reactant_ids_bytes);  
        cudaMalloc(&d_number_of_products_, number_of_products_bytes);
        cudaMalloc(&d_product_ids_, product_ids_bytes);  
        cudaMalloc(&d_yields_, yields_bytes); 

        //copy data from host memory to device memory    
        cudaMemcpy(d_rate_constants, rate_constants_data, rate_constants_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_state_variables, state_variables_data, state_forcing_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_forcing, forcing_data, state_forcing_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_number_of_reactants_, number_of_reactants, number_of_reactants_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_reactant_ids_, reactant_ids, reactant_ids_bytes,cudaMemcpyHostToDevice); 
        cudaMemcpy(d_number_of_products_, number_of_products, number_of_products_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_product_ids_, product_ids, product_ids_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_yields_, yields, yields_bytes, cudaMemcpyHostToDevice); 

        //total thread count == number of grid cells
        int block_size = 32; 
        int num_block = (ngrids + block_size - 1)/block_size; 
        
        //kernel function call
        AddForcingTerms_kernel<<<num_block, block_size>>>(
            d_rate_constants, 
            d_state_variables, 
            d_forcing, 
            ngrids, nrxns, nspecs,
            d_number_of_reactants_, 
            d_reactant_ids_, 
            d_number_of_products_, 
            d_product_ids_, 
            d_yields_);
        cudaDeviceSynchronize(); 
        
        // copy data from device memory to host memory
        cudaMemcpy(forcing_data, d_forcing, state_forcing_bytes, cudaMemcpyDeviceToHost);
   
        // clean up	
        cudaFree(d_rate_constants); 
        cudaFree(d_state_variables); 
        cudaFree(d_forcing);
        cudaFree(d_number_of_reactants_);
        cudaFree(d_reactant_ids_);
        cudaFree(d_number_of_products_);
        cudaFree(d_product_ids_);
        cudaFree(d_yields_);
    } // end of AddForcingTerms_kernelSetup
    } // namespace cuda 
}     // namespace micm
