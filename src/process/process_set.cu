#include <iostream>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <iostream>


namespace micm {
    namespace cuda {
    //one thread per reaction
    //passing all device pointers 

    __global__ void AddForcingTerms_kernel(
        double* rate_constants, 
        double* state_variables, 
        double* forcing, 
        int n_grids, 
        int n_reactions, 
        int n_species,
        size_t* number_of_reactants_, 
        size_t* accumulated_n_reactants, 
        size_t* reactant_ids_, 
        size_t* number_of_products_, 
        size_t* accumulated_n_products, 
        size_t* product_ids_, 
        double* yields_)
    {
    //define thread index 
    //one thread per reaction
    int tid = blockIdx.x * blockDim.x + threadIdx.x; 
    int rate_constants_size = n_grids * n_reactions; 
    
    if (tid < rate_constants_size){
        double rate = rate_constants[tid];
        int grid_index = tid % n_grids; 
        int reaction_index = (tid - grid_index)/n_grids; 

        int reactant_num = number_of_reactants_[reaction_index]; //number of reactants of the reaction
        int product_num = number_of_products_[reaction_index]; //number of products of the reaction 
        
        int initial_reactant_ids_index = accumulated_n_reactants[reaction_index];
        int initial_product_ids_index = accumulated_n_products[reaction_index];
        int initial_yields_index = accumulated_n_products[reaction_index]; 
        
        for (int i_reactant = 0; i_reactant < reactant_num; i_reactant++){
            rate *= state_variables[reactant_ids_[initial_reactant_ids_index + i_reactant] * n_grids + grid_index];  
        }
        
        for (int i_reactant = 0; i_reactant < reactant_num; i_reactant++){
            double rate_subtration = 0 - rate; 
            atomicAdd(&forcing[reactant_ids_[initial_reactant_ids_index + i_reactant] * n_grids + grid_index], rate_subtration); 
        }
        for (int i_product = 0; i_product < product_num; i_product++){
            atomicAdd(&forcing[product_ids_[initial_product_ids_index + i_product] * n_grids + grid_index] , yields_[initial_yields_index + i_product] * rate); 
        
        } //looping number of product times
    } // checking valid tid value
  } //AddForcingTerms_kernel function
    
    void AddForcingTerms_kernelSetup(
        const size_t* number_of_reactants,
        int number_of_reactants_size,
        const size_t* reactant_ids, 
        int reactant_ids_size,
        const size_t* number_of_products, 
        int number_of_products_size,
        const size_t* product_ids,
        int product_ids_size,
        const double* yields,
        int yields_size,
        const Matrix<double>& rate_constants, 
        const Matrix<double>& state_variables, 
        Matrix<double>& forcing)
    {
        //data of matrices
        int n_grids = rate_constants[0].size();
        int n_reactions = rate_constants.size(); 
        int n_species = state_variables.size(); 

        const double* rate_constants_data = rate_constants.AsVector().data(); 
        const double* state_variables_data = state_variables.AsVector().data();
        double* forcing_data = forcing.AsVector().data(); 
       
       
        int accumulated_n_reactants_bytes = sizeof(size_t) * (number_of_reactants_size); 
        size_t* accumulated_n_reactants = (size_t*)malloc(accumulated_n_reactants_bytes); 
        accumulated_n_reactants[0] = 0; 
        
        for (int i = 0; i < number_of_reactants_size - 1; i++){
            int sum = accumulated_n_reactants[i] + number_of_reactants[i]; 
            accumulated_n_reactants[i+1] = sum; 
        } 

        int accumulated_n_products_bytes = sizeof(size_t) * (number_of_products_size); 
        size_t* accumulated_n_products = (size_t*)malloc(accumulated_n_products_bytes); 
        accumulated_n_products[0] = 0;  
        
        for (int k = 0; k < number_of_products_size - 1; k++){
            int sum = accumulated_n_products[k] + number_of_products[k]; 
            accumulated_n_products[k+1] = sum; 
        }

        // device pointer to vectors
        double* d_rate_constants; 
        double* d_state_variables; 
        double* d_forcing; 
        size_t* d_number_of_reactants_; 
        size_t* d_accumulated_n_reactants; 
        size_t* d_reactant_ids_; 
        size_t* d_number_of_products_; 
        size_t* d_accumulated_n_products; 
        size_t* d_product_ids_; 
        double* d_yields_; 
       
    
        //allocate device memory
        size_t rate_constants_bytes = sizeof(double) * (n_grids * n_reactions); 
        size_t state_forcing_bytes = sizeof(double) * (n_grids * n_species); 
        size_t number_of_reactants_bytes = sizeof(size_t) * number_of_reactants_size;
        size_t reactant_ids_bytes = sizeof(size_t) * reactant_ids_size; 
        size_t number_of_products_bytes = sizeof(size_t) * number_of_products_size; 
        size_t product_ids_bytes = sizeof(size_t) * product_ids_size;
        size_t yields_bytes = sizeof(double) * yields_size;
        
        cudaMalloc(&d_rate_constants, rate_constants_bytes); 
        cudaMalloc(&d_state_variables, state_forcing_bytes); 
        cudaMalloc(&d_forcing, state_forcing_bytes); 
        cudaMalloc(&d_number_of_reactants_, number_of_reactants_bytes);
        cudaMalloc(&d_accumulated_n_reactants, accumulated_n_reactants_bytes); 
        cudaMalloc(&d_reactant_ids_, reactant_ids_bytes);  
        cudaMalloc(&d_number_of_products_, number_of_products_bytes);
        cudaMalloc(&d_accumulated_n_products, accumulated_n_products_bytes); 
        cudaMalloc(&d_product_ids_, product_ids_bytes);  
        cudaMalloc(&d_yields_, yields_bytes); 

        //copy data from host memory to device memory    
        cudaMemcpy(d_rate_constants, rate_constants_data, rate_constants_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_state_variables, state_variables_data, state_forcing_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_forcing, forcing_data, state_forcing_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_number_of_reactants_, number_of_reactants, number_of_reactants_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_accumulated_n_reactants, accumulated_n_reactants, accumulated_n_reactants_bytes,cudaMemcpyHostToDevice); 
        cudaMemcpy(d_reactant_ids_, reactant_ids, reactant_ids_bytes,cudaMemcpyHostToDevice); 
        cudaMemcpy(d_number_of_products_, number_of_products, number_of_products_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_accumulated_n_products, accumulated_n_products, accumulated_n_products_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_product_ids_, product_ids, product_ids_bytes, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_yields_, yields, yields_bytes, cudaMemcpyHostToDevice); 

        //total thread count == rate_constants matrix size
        int N = n_grids * n_reactions; 
        int block_size = 320; 
        int num_block = (N + block_size -1)/block_size; 
        
        //kernel function call
        AddForcingTerms_kernel<<<num_block, block_size>>>( 
            d_rate_constants, 
            d_state_variables, 
            d_forcing, 
            n_grids, 
            n_reactions, 
            n_species, 
            d_number_of_reactants_, 
            d_accumulated_n_reactants, 
            d_reactant_ids_, 
            d_number_of_products_, 
            d_accumulated_n_products, 
            d_product_ids_, 
            d_yields_);
        cudaDeviceSynchronize(); 
        
        cudaMemcpy(forcing_data, d_forcing, state_forcing_bytes, cudaMemcpyDeviceToHost);
        
        //clean up
        cudaFree(d_rate_constants); 
        cudaFree(d_state_variables); 
        cudaFree(d_forcing);
        cudaFree(d_number_of_reactants_);
        cudaFree(d_accumulated_n_reactants);
        cudaFree(d_reactant_ids_);
        cudaFree(d_number_of_products_);
        cudaFree(d_accumulated_n_products);
        cudaFree(d_product_ids_);
        cudaFree(d_yields_ );
        free(accumulated_n_reactants); 
        free(accumulated_n_products); 
    
        
     } //kernel function setup
    }//namespace cuda 
}//namespace micm
