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


__global__ void AddJacobianTerms_kernel(
  double* rate_constants,
  double* state_variables,
  double* jacobian,
  int n_grids,
  int n_reactions,
  int n_species,
  size_t*  number_of_reactants,
  size_t* reactant_ids,
  size_t* number_of_products,
  size_t* product_ids,
  double* yields,
  size_t* jacobian_flat_ids,
  size_t rows_ids_size){
  
  int tid = blockIdx.x * blockDim.x + threadIdx.x; 
  size_t react_id_offset, prod_id_offset, yield_offset, jacobian_id_idx, initial_jacobian_idx; 
   
    //paralle grids -> one thread per grid cell
    //how to get jacobian initial index for each grid
    if (tid < n_grids){
      printf("running in kernel\n"); 
      react_id_offset = 0, prod_id_offset = 0, yield_offset = 0, jacobian_id_idx = -1; 
      size_t* jacobian_flat_ids_ptr = jacobian_flat_ids; 
      initial_jacobian_idx = tid * rows_ids_size; 
      for (int i_rxn = 0; i_rxn < n_reactions; i_rxn++){
          printf ("reaction index %d\n",i_rxn); 
          for(int i_ind = 0; i_ind < number_of_reactants[i_rxn]; i_ind++){
             double d_rate_d_int = rate_constants[i_rxn * n_grids + tid]; 
             for (int i_react = 0; i_react < number_of_reactants[i_rxn]; i_react++){
              if (i_react != i_ind){
                d_rate_d_int *= state_variables[reactant_ids[react_id_offset + i_react] * n_grids + tid];
              }
             }//second inner loop
       }//first inner loop
      }//loop over num_reactions
    }//check for valid tid 
  }// end of AddJacobianTerms_kernel
    

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
    

    void AddJacobianTerms_kernelSetup(
        const double* rate_constants,
        const double* state_variables,
        double* jacobian,
        size_t n_grids,
        size_t n_reactions,
        size_t n_species,
        size_t jacobian_size,
        size_t row_ids_size,
        const size_t* number_of_reactants, 
        const size_t* reactant_ids, 
        size_t reactant_ids_size, 
        const size_t* number_of_products, 
        const size_t* product_ids,
        size_t product_ids_size, 
        const double* yields,
        size_t yields_size,
        size_t* jacobian_flat_ids,
        size_t jacobian_flat_ids_size){
        std::cout << "grid size: "<< n_grids<<std::endl; 
        std::cout << "reaction size: "<< n_reactions<<std::endl; 
        std::cout << "species size: "<<n_species<<std::endl; 
        //device pointer
        double* d_rate_constants; 
        double* d_state_variables; 
        double* d_jacobian; 
        size_t* d_number_of_reactants;    
        size_t* d_reactant_ids; 
        size_t* d_number_of_products; 
        size_t* d_product_ids; 
        double* d_yields; 
        size_t* d_jacobian_flat_ids; 

        //allocate device memory 
        cudaMalloc(&d_rate_constants, sizeof(double)* n_grids * n_reactions); 
        cudaMalloc(&d_state_variables, sizeof(double)* n_grids * n_species); 
        cudaMalloc(&d_jacobian, sizeof(double)* jacobian_size);
        cudaMalloc(&d_number_of_reactants, sizeof(size_t)* n_reactions); 
        cudaMalloc(&d_reactant_ids, sizeof(size_t) * reactant_ids_size);
        cudaMalloc(&d_number_of_products, sizeof(size_t)* n_reactions);
        cudaMalloc(&d_product_ids, sizeof(size_t) * product_ids_size);  
        cudaMalloc(&yields, sizeof(double) * yields_size); 
        cudaMalloc(&jacobian_flat_ids, sizeof(size_t)* jacobian_flat_ids_size); 


        //transfer data from host to device 
        cudaMemcpy(d_rate_constants, rate_constants, sizeof(double)* n_grids * n_reactions,cudaMemcpyHostToDevice); 
        cudaMemcpy(d_state_variables, state_variables, sizeof(double)* n_grids * n_species, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_jacobian, jacobian, sizeof(double)* jacobian_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_number_of_reactants, number_of_reactants, sizeof(size_t)* n_reactions,cudaMemcpyHostToDevice); 
        cudaMemcpy(d_reactant_ids, reactant_ids, sizeof(size_t) * reactant_ids_size, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_number_of_products, number_of_products, sizeof(size_t)* n_reactions, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_product_ids, product_ids, sizeof(size_t)* product_ids_size, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_yields, yields, sizeof(double) * yields_size, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_jacobian_flat_ids, jacobian_flat_ids, sizeof(size_t)* jacobian_flat_ids_size, cudaMemcpyHostToDevice); 

        //total thread count == n_grids 
        int block_size = 32; 
        int num_blocks = (n_grids + block_size -1)/block_size; 
        //kernel function call
        AddJacobianTerms_kernel<<<num_blocks, block_size>>>(
          d_rate_constants,
          d_state_variables,
          d_jacobian,
          n_grids,
          n_reactions,
          n_species,
          d_number_of_reactants, 
          d_reactant_ids, 
          d_number_of_products, 
          d_product_ids,
          d_yields,
          d_jacobian_flat_ids,
          row_ids_size
        );
          cudaDeviceSynchronize(); 
          cudaMemcpy(jacobian, d_jacobian, sizeof(double)* jacobian_size, cudaMemcpyDeviceToHost);
        
        // clean up	
        cudaFree(d_rate_constants); 
        cudaFree(d_state_variables); 
        cudaFree(d_jacobian);
        cudaFree(d_number_of_reactants);
        cudaFree(d_reactant_ids);
        cudaFree(d_number_of_products);
        cudaFree(d_product_ids);
        cudaFree(d_yields);
        cudaFree(d_jacobian_flat_ids); 
    }//end of AddJacobianTerms_kernelSetup function 
  } // namespace cuda 
}     // namespace micm
