#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <iostream>

namespace micm {
    namespace cuda {
    //one thread per reaction
    //passing all device pointers 
    __global__ void AddForcingTerms_kernel(
        double* rate_array,
        double* state_variable, 
        double* rate_constants, 
        double* state_variables, 
        double* forcing, 
        int matrix_rows, 
        int rate_constants_columns, 
        int state_forcing_columns,
        size_t* number_of_reactants_, 
        size_t* accumulated_n_reactants, 
        size_t* reactant_ids_, 
        size_t* number_of_products_, 
        size_t* accumulated_n_products, 
        size_t* product_ids_, 
        size_t* yields_)
    {
    //define thread index 
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x; 
    int rate_constants_size = matrix_rows * rate_constants_columns; 
    
    if (tid < rate_constants_size){
        int rate = rate_constants[tid];
        int rate_constants_col_index = tid % rate_constants_columns; 
        int row_index = (tid - rate_constants_col_index)/rate_constants_columns;
    
        int reactant_num = number_of_reactants_[rate_constants_col_index]; //number of reactants of the reaction
        int product_num = number_of_products_[rate_constants_col_index]; //number of products of the reaction 
        
        int initial_reactant_ids_index = accumulated_n_reactants[rate_constants_col_index];
        int initial_product_ids_index = accumulated_n_products[rate_constants_col_index];
        int initial_yields_index = accumulated_n_products[rate_constants_col_index]; 
        
       
        for (int i_reactant = 0; i_reactant < reactant_num; i_reactant++){
            int reactant_ids_index = initial_reactant_ids_index + i_reactant; 
            int state_forcing_col_index = reactant_ids_[reactant_ids_index]; 
            //debugging
            if (tid == 0){
                state_variable[tid + i_reactant] = state_variables[row_index * state_forcing_columns + state_forcing_col_index];  
            }
            rate *= state_variables[row_index * state_forcing_columns + state_forcing_col_index];  
        }
        //debugging 
        rate_array[tid] = rate; 
        
        for (int i_reactant = 0; i_reactant < reactant_num; i_reactant++){
            int reactant_ids_index = initial_reactant_ids_index + i_reactant; 
            int state_forcing_col_index = reactant_ids_[reactant_ids_index]; 
            forcing[row_index * state_forcing_columns + state_forcing_col_index] -=rate; 
        }
        
        for (int i_product = 0; i_product < product_num; i_product++){
            int yields_index = initial_yields_index + i_product; 
            int product_ids_index  = initial_product_ids_index + i_product; 
            int forcing_col_index = product_ids_[product_ids_index]; 
            forcing[row_index * state_forcing_columns + forcing_col_index] += yields_[yields_index] * rate; 
        } 
        
    }
  }
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
        int matrix_rows = rate_constants.size(); 
        int rate_constants_columns = rate_constants[0].size(); 
        int state_forcing_columns = state_variables[0].size();

        const double* rate_constants_data = rate_constants.AsVector().data(); 
        const double* state_variables_data = state_variables.AsVector().data();
        double* forcing_data = forcing.AsVector().data(); 
       
       //this vector provides initial index to reactant_ids_ to get reactant id of every reaction 
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
        size_t* d_yields_; 
       
        //debugging
        int rate_array_size = matrix_rows * rate_constants_columns; 
        double* d_rate_array;
        double* rate_array; 
        cudaMalloc(&d_rate_array, sizeof(double) * rate_array_size); 
        rate_array = (double*)malloc(sizeof(double) * rate_array_size);
        double* d_state_variable; 
        double* state_variable;
        cudaMalloc(&d_state_variable, 2); 
        rate_array = (double*)malloc(sizeof(double) * 2);
        
        //allocate device memory
        size_t rate_constants_bytes = sizeof(double) * (matrix_rows * rate_constants_columns); 
        size_t state_forcing_bytes = sizeof(double) * (matrix_rows * state_forcing_columns); 
        size_t number_of_reactants_bytes = sizeof(size_t) * number_of_reactants_size;
        size_t reactant_ids_bytes = sizeof(size_t) * reactant_ids_size; 
        size_t number_of_products_bytes = sizeof(size_t) * number_of_products_size; 
        size_t product_ids_bytes = sizeof(size_t) * product_ids_size;
        size_t yields_bytes = sizeof(size_t) * yields_size;
        
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
        int N = matrix_rows * rate_constants_columns; 
        int block_size = 320; 
        int num_block = (N + block_size -1)/block_size; 
        //kernel function call
    
        AddForcingTerms_kernel<<<num_block, block_size>>>(
            d_rate_array,
            d_state_variable,
            d_rate_constants, 
            d_state_variables, 
            d_forcing, 
            matrix_rows, 
            rate_constants_columns, 
            state_forcing_columns, 
            d_number_of_reactants_, 
            d_accumulated_n_reactants, 
            d_reactant_ids_, 
            d_number_of_products_, 
            d_accumulated_n_products, 
            d_product_ids_, 
            d_yields_);
        cudaDeviceSynchronize(); 
        cudaMemcpy(forcing_data, d_forcing, state_forcing_bytes, cudaMemcpyDeviceToHost);
        
        //debugging 
        cudaMemcpy(rate_array, d_rate_array, sizeof(double)*rate_array_size, cudaMemcpyDeviceToHost );    
        std::cout << "this is rate_array: "<< std::endl; 
        for (int k = 0; k < rate_array_size; k++){
            std::cout << rate_array[k]<<std::endl; 
        }
        cudaMemcpy(state_variable, d_state_variable, sizeof(double)*2, cudaMemcpyDeviceToHost);   
        for (int k = 0; k < 2; k++){
            std::cout << state_variable[k]<<std::endl; 
        }
       
       
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
        }
    }//namespace cuda 
}//namespace micm