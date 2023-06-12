#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <vector>


 void AddForcingTerms_kernelSetup(
      const Matrix<double>& rate_constants,
      const Matrix<double>& state_variables,
      Matrix<double>& forcing, vector<std::size_t> number_of_reactants_, 
      vector<std::size_t> reactant_ids_, vector<std::size_t> number_of_products_,
      vector<std::size_t> product_ids_, vector<std::size_t> yields_) const
    {
        vector <std::size_t> accumulated_n_reactants; 
        accumulated_n_reactants.push_back(0); 
        for (int i = 0; i < number_of_reactants_.size(); i++){
            int sum = accumulated_n_reactants[i] + number_of_reactants_[i]; 
            accumulated_n_reactants.push_back(sum); 
        }
        
        vector<std::size_t>accumulated_n_products;
        accumulated_n_products.push_back(0); 
        
        for (int i = 0; i < number_of_products_.size(); i++){
            int sum = accumulated_n_reactants[i] + number_of_products_[i]; 
            accumlated_n_products.push_back[sum];     
        }

        // device pointer to vectors
        std::vector<double>* d_rate_constants; 
        std::vector<double>* d_state_variables; 
        std::vector<double>* d_forcing; 
        std::vector<std::size_t>* d_number_of_reactants_; 
        std::vector<std::size_t>* d_accumulated_n_reactants; 
        std::vector<std::size_t>* d_reactant_ids_; 
        std::vector<std::size_t>* d_number_of_products_; 
        std::vector<std::size_t>* d_accumulated_n_products; 
        std::vector<std::size_t>* d_product_ids_; 
        std::vector<std::size_t>* d_yields_; 
        


        //allocate device memory
        cudaMalloc(&d_rate_constants, (sizeof(std::double)* rate_constants.data.size())); 
        cudaMalloc(&d_state_variables, (sizeof(std::double)* state_variables.data.size())); 
        cudaMalloc(&d_forcing, (sizeof(std::double)* forcing.data.size())); 
        cudaMalloc(&d_number_of_reactants_, (sizeof(std::size_t)* number_of_reactants_.size()));
        cudaMalloc(&d_accumulated_n_reactants,(sizeof(std::size_t)* accumulated_n_reactants.size())); 
        cudaMalloc(&d_reactant_ids_, (sizeof(std::size_t)* reactant_ids_.size()));  
        cudaMalloc(&d_number_of_products_, (sizeof(std::size_t)* number_of_products_.size()));
        cudaMalloc(&d_accumulated_n_products, (sizeof(std::size_t) * accumulated_n_products.size())); 
        cudaMalloc(&d_product_ids_, (sizeof(std::size_t)* product_ids_.size()));  
        cudaMalloc(&d_yields, (sizeof(std::size_t) * yields_.size())); 
     
       //copy data from host to device memory 
        cudaMemcpy(d_rate_constants, &rate_constants.data, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_state_variables, &state_variables.data, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_forcing, &forcing.data, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_number_of_reactants_, &number_of_reactants_, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_accumulated_n_reactants, &accumulated_n_reactants, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_reactant_ids_, &reactant_ids_, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_number_of_products_, &number_of_products_, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_accumulated_n_products, &accumlated_n_products, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_product_ids_, &product_ids_, cudaMemcpyHostToDevice); 
        cudaMemcpy(d_yields, &yields_, cudaMemcpyHostToDevice); 

        //total thread count == rate_constants matrix size?
        int threads_count = rate_constants.x_dim * rate_constants.y_dim;
        //block size 
        int threadsPerBlock = 128; //32 threads per warp * 4 warps
        //grid size 
        int blocks_count = (int)ceil(threads_count/threadsPerBlock); 
    
        int matrix_rows = rate_constants.x_dim; 
        int rate_constants_columns = rate_constants.y_dim; 
        int state_forcing_columns = state_variables.y_dim;

        //kernel function call
        AddForcingTerms_kernel(d_rate_constants, d_state_variables, 
        d_forcing, matrix_rows, rate_constants_columns, state_forcing_columns, 
        number_of_reactants_, accumulated_n_reactants, reactant_ids_, 
        number_of_products_, accumulated_n_products, product_ids_, 
        yields_);
        cudaDeviceSynchronize(); 
        
        cudaMemcpy(d_forcing, &forcing.data, cudaMemcpyDeviceToHost);

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
    
    }

  
//one thread per reaction in atompheric model 
__global__ void AddForcingTerms_kernel(vector<std:: double>* rate_constants, vector<std::double> state_variables, 
    vector<std::double> forcing, int matrix_rows, int rate_constants_columns, int state_forcing_columns, 
    vector<std::size_t> number_of_reactants_, vector<std::size_t> accumulated_n_reactants, vector<std::size_t> reactant_ids_, 
    vector<std::size_t> number_of_products_,vector<std::size_t> accumulated_n_products, vector<std::size_t> product_ids_, 
    vector<std::size_t> yields_){

    //define thread index 
    int tid = blockIdx.x + blockDim.x + threadIdx.x; 
   
    if (tid < rate_constants.size()){
        int rate = rate_constants[tid]; // rate of a specific reaction in a specific gridcell 
        int row_index = tid % rate_constants_columns; 
        int reactant_num = number_of_reactants_[tid % rate_constants_columns]; //number of reactants of the reaction
        int product_num = number_of_products_[tid % rate_constants_columns]; //number of products of the reaction 
        int initial_reactant_ids_index = accumulated_n_reactants[tid % rate_constants_columns];
        int initial_product_ids_index = accumulated_n_products[tid % rate_constants_columns];
        int initial_yields_index = accumulated_n_products[tid % rate_constants_columns]; 
        
        //access index at reactant_ids based on number_of_reactant_
        for (int i_reactant = 0; i_reactant < reactant_num; i_reactant++){
            int reactant_ids_index = i_reactant + initial_reactant_ids_index; 
            int state_forcing_col_index = reactant_ids_[reactant_ids_index]; 
            //how to match thread idx to state_variable index 
            //but we need to consider the row of state_variable 
            rate *= state_variables[row_index * state_forcing_columns + state_forcing_col_index]; 
            forcing[row_index * state_forcing_columns + state_forcing_col_index] -= rate; 
        }
        for (int i_product = 0; i_product < product_num; i_product++){
            int yields_index = initial_yields_index + i_product; 
            int product_ids_index  = initial_product_ids_index + i_product; 
            int forcing_col_index = product_ids_[product_ids_index]; 
            forcing[row_index * state_forcing_columns + forcing_col_index] += yields_[yields_index] * rate; 
        }   
    }
}

