#include <iostream>

namespace micm
{
  namespace cuda
  {
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
      // define thread index
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      size_t react_id_offset, prod_id_offset, yield_offset;

      if (tid < ngrids)
      {
        react_id_offset = 0;
        prod_id_offset = 0;
        yield_offset = 0;
        for (std::size_t i_rxn = 0; i_rxn < nrxns; ++i_rxn)
        {
          double rate = rate_constants[i_rxn * ngrids + tid];
          for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
            rate *= state_variables[reactant_ids_[react_id_offset + i_react] * ngrids + tid];
          for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
          {
            forcing[reactant_ids_[react_id_offset + i_react] * ngrids + tid] -= rate;
          }
          for (std::size_t i_prod = 0; i_prod < number_of_products_[i_rxn]; ++i_prod)
          {
            size_t index = product_ids_[prod_id_offset + i_prod] * ngrids + tid;
            forcing[index] += yields_[yield_offset + i_prod] * rate;
          }
          react_id_offset += number_of_reactants_[i_rxn];
          prod_id_offset += number_of_products_[i_rxn];
          yield_offset += number_of_products_[i_rxn];
        }  // for loop over number of reactions
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
  size_t* acc_n_reactants,
  size_t* reactant_ids,
  size_t* number_of_products,
  size_t* acc_n_products,
  size_t* product_ids,
  double* yields,
  size_t* jacobian_flat_ids,
  size_t* acc_n_jacobian_flat_ids,
  size_t row_ids_size){
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x; 
    int rate_constants_size = n_grids*n_reactions; 
    
    if (tid < rate_constants_size){
    
    double d_rate_d_ind = rate_constants[tid]; 
    int reaction_idx = tid % n_reactions; 
    int grid_idx = (tid - reaction_idx)/n_reactions; 
    int num_reactants = number_of_reactants[reaction_idx]; 
    int num_products = number_of_products[reaction_idx];
    int initial_reactant_ids_idx = acc_n_reactants[reaction_idx]; 
    int initial_yields_idx = acc_n_products[reaction_idx]; 
    int initial_jacobian_flat_ids_idx = acc_n_jacobian_flat_ids[reaction_idx]; 
    int acc_jacobian_flat_ids_idx = 0; 
    int initial_jacobian_idx = grid_idx * row_ids_size; 
    printf("tid: %d\n, num_reactants: %d\n", tid, num_reactants); 

    //loop over over num_reactants of every reaction
    for (int i_ind = 0; i_ind < num_reactants; i_ind++){
        
        for (int i_react = 0; i_react < num_reactants; i_react){
          if (i_ind != i_react){
            d_rate_d_ind *= state_variables[grid_idx * n_species + reactant_ids[initial_reactant_ids_idx + i_react]];
          }
        }
        for(int i_dep = 0; i_dep < num_reactants;i_dep++){
          int jacobian_idx = initial_jacobian_idx + jacobian_flat_ids[initial_jacobian_flat_ids_idx+i_dep];
          acc_jacobian_flat_ids_idx = initial_jacobian_flat_ids_idx+i_dep;
          jacobian[jacobian_idx] -= d_rate_d_ind; 
        }
        for (int i_dep = 0; i_dep < num_products; i_dep++){
          int jacobian_idx = initial_jacobian_idx + jacobian_flat_ids[acc_jacobian_flat_ids_idx + i_dep]; 
          jacobian[jacobian_idx] += yields[initial_yields_idx + i_dep] * d_rate_d_ind;
        }
      }
    }
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
        int yields_size)
    {
      // device pointer to vectorss
      double* d_rate_constants;
      double* d_state_variables;
      double* d_forcing;
      double* d_yields_;
      size_t* d_number_of_reactants_;
      size_t* d_reactant_ids_;
      size_t* d_number_of_products_;
      size_t* d_product_ids_;

      // allocate device memory
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

      // copy data from host memory to device memory
      cudaMemcpy(d_rate_constants, rate_constants_data, rate_constants_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(d_state_variables, state_variables_data, state_forcing_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(d_forcing, forcing_data, state_forcing_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(d_number_of_reactants_, number_of_reactants, number_of_reactants_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(d_reactant_ids_, reactant_ids, reactant_ids_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(d_number_of_products_, number_of_products, number_of_products_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(d_product_ids_, product_ids, product_ids_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(d_yields_, yields, yields_bytes, cudaMemcpyHostToDevice);

      // total thread count == number of grid cells
      int block_size = 32;
      int num_block = (ngrids + block_size - 1) / block_size;

      // kernel function call
      AddForcingTerms_kernel<<<num_block, block_size>>>(
          d_rate_constants,
          d_state_variables,
          d_forcing,
          ngrids,
          nrxns,
          nspecs,
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
    }  // end of AddForcingTerms_kernelSetup
  }    // namespace cuda
}  // namespace micm
