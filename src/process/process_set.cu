#include <iostream>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  // one thread per reaction
  // passing all device pointers
  __global__ void AddForcingTerms_kernel(
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
      double* yields_)
  {
    // define thread index
    int tid = blockIdx.x + blockDim.x + threadIdx.x;
    int rate_reactants_size = matrix_rows * rate_constants_columns;
    if (tid < rate_reactants_size)
    {
      double rate = rate_constants[tid];  // rate of a specific reaction in a specific gridcell
      int row_index = tid % rate_constants_columns;
      int reactant_num = number_of_reactants_[tid % rate_constants_columns];  // number of reactants of the reaction
      int product_num = number_of_products_[tid % rate_constants_columns];    // number of products of the reaction
      int initial_reactant_ids_index = accumulated_n_reactants[tid % rate_constants_columns];
      int initial_product_ids_index = accumulated_n_products[tid % rate_constants_columns];
      int initial_yields_index = accumulated_n_products[tid % rate_constants_columns];

      // access index at reactant_ids based on number_of_reactant_
      for (int i_reactant = 0; i_reactant < reactant_num; i_reactant++)
      {
        int reactant_ids_index = i_reactant + initial_reactant_ids_index;
        int state_forcing_col_index = reactant_ids_[reactant_ids_index];
        // how to match thread idx to state_variable index
        // but we need to consider the row of state_variable
        rate *= state_variables[row_index * state_forcing_columns + state_forcing_col_index];
      }

      for (int i_reactant = 0; i_reactant < reactant_num; i_reactant++)
      {
        int reactant_ids_index = i_reactant + initial_reactant_ids_index;
        int state_forcing_col_index = reactant_ids_[reactant_ids_index];

        forcing[row_index * state_forcing_columns + state_forcing_col_index] -= rate;
      }

      for (int i_product = 0; i_product < product_num; i_product++)
      {
        int yields_index = initial_yields_index + i_product;
        int product_ids_index = initial_product_ids_index + i_product;
        int forcing_col_index = product_ids_[product_ids_index];
        forcing[row_index * state_forcing_columns + forcing_col_index] += yields_[yields_index] * rate;
      }
    }
  }

  void AddForcingTerms_kernelSetup(
      const Matrix<double>& rate_constants,
      const Matrix<double>& state_variables,
      Matrix<double>& forcing,
      const std::vector<std::size_t>& number_of_reactants_,
      const std::vector<std::size_t>& reactant_ids_,
      const std::vector<std::size_t>& number_of_products_,
      const std::vector<std::size_t>& product_ids_,
      const std::vector<double>& yields_,
      const std::vector<std::size_t>& jacobian_flat_ids_)
  {
    int matrix_rows = rate_constants.size();
    int rate_constants_columns = rate_constants[0].size();
    int state_forcing_columns = state_variables[0].size();
    // access class member vectors
    const size_t* number_of_reactants = number_of_reactants_.data();
    int number_of_reactants_size = number_of_reactants_.size();
    const size_t* reactant_ids = reactant_ids_.data();
    int reactant_ids_size = reactant_ids_.size();
    const size_t* number_of_products = number_of_products_.data();
    int number_of_products_size = number_of_products_.size();
    const size_t* product_ids = product_ids_.data();
    int product_ids_size = product_ids_.size();
    const double* yields = yields_.data();
    int yields_size = yields_.size();

    // allocate memory for host c pointers
    int accumulated_n_reactants_bytes = sizeof(size_t) * (number_of_reactants_size + 1);
    size_t* accumulated_n_reactants = (size_t*)malloc(accumulated_n_reactants_bytes);
    accumulated_n_reactants[0] = 0;
    for (int i = 0; i < number_of_reactants_size; i++)
    {
      int sum = accumulated_n_reactants[i] + number_of_reactants[i];
      accumulated_n_reactants[i + 1] = sum;
    }

    int accumulated_n_products_bytes = sizeof(size_t) * (number_of_products_size + 1);
    size_t* accumulated_n_products = (size_t*)malloc(accumulated_n_products_bytes);
    accumulated_n_products[0] = 0;
    for (int i = 0; i < number_of_products_size; i++)
    {
      int sum = accumulated_n_products[i] + number_of_products[i];
      accumulated_n_products[i + 1] = sum;
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

    // allocate device memory
    size_t rate_constants_bytes = sizeof(double) * (matrix_rows * rate_constants_columns);
    size_t state_variables_bytes = sizeof(double) * (matrix_rows * state_forcing_columns);
    size_t forcing_bytes = sizeof(double) * (matrix_rows * state_forcing_columns);
    size_t number_of_reactants_bytes = sizeof(size_t) * number_of_reactants_size;
    size_t reactant_ids_bytes = sizeof(size_t) * reactant_ids_size;
    size_t number_of_products_bytes = sizeof(size_t) * number_of_products_size;
    size_t product_ids_bytes = sizeof(size_t) * product_ids_size;
    size_t yields_bytes = sizeof(double) * yields_size;

    cudaMalloc(&d_rate_constants, rate_constants_bytes);
    cudaMalloc(&d_state_variables, state_variables_bytes);
    cudaMalloc(&d_forcing, forcing_bytes);
    cudaMalloc(&d_number_of_reactants_, number_of_reactants_bytes);
    cudaMalloc(&d_accumulated_n_reactants, accumulated_n_reactants_bytes);
    cudaMalloc(&d_reactant_ids_, reactant_ids_bytes);
    cudaMalloc(&d_number_of_products_, number_of_products_bytes);
    cudaMalloc(&d_accumulated_n_products, accumulated_n_products_bytes);
    cudaMalloc(&d_product_ids_, product_ids_bytes);
    cudaMalloc(&d_yields_, yields_bytes);

    // copy data from host to device memory
    const double* rate_constants_data = rate_constants.AsVector().data();
    const double* state_variables_data = state_variables.AsVector().data();
    double* forcing_data = forcing.AsVector().data();

    cudaMemcpy(d_rate_constants, rate_constants_data, rate_constants_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_state_variables, state_variables_data, state_variables_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_forcing, forcing_data, forcing_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_number_of_reactants_, number_of_reactants, number_of_reactants_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_accumulated_n_reactants, accumulated_n_reactants, accumulated_n_reactants_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_reactant_ids_, reactant_ids, reactant_ids_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_number_of_products_, number_of_products, number_of_products_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_accumulated_n_products, accumulated_n_products, accumulated_n_products_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_product_ids_, product_ids, product_ids_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_yields_, yields, yields_bytes, cudaMemcpyHostToDevice);

    // total thread count == rate_constants matrix size?
    int threads_count = matrix_rows * rate_constants_columns;
    // block size
    int threadsPerBlock = 128;  // 32 threads per warp * 4 warps
    // grid size
    int blocks_count = (int)ceil(threads_count / threadsPerBlock);

    // kernel function call
    AddForcingTerms_kernel<<<blocks_count, threadsPerBlock>>>(
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
    cudaMemcpy(d_forcing, forcing_data, forcing_bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_rate_constants);
    cudaFree(d_state_variables);
    cudaFree(d_forcing);
    cudaFree(d_number_of_reactants_);
    cudaFree(d_accumulated_n_reactants);
    cudaFree(d_reactant_ids_);
    cudaFree(d_number_of_products_);
    cudaFree(d_accumulated_n_products);
    cudaFree(d_product_ids_);
    cudaFree(d_yields_);

    free(accumulated_n_reactants);
    free(accumulated_n_products);
  }
}  // namespace micm