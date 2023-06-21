#pragma once


#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <vector>

namespace micm {
    namespace cuda {

    void AddForcingTerms_kernelSetup(
        size_t* number_of_reactants,
        int number_of_reactants_size,
        size_t* reactant_ids, 
        int reactant_ids_size,
        size_t* number_of_products, 
        int number_of_products_size,
        size_t* product_ids,
        int product_ids_size,
        double* yields,
        int yields_size,
        const Matrix<double>& rate_constants, 
        const Matrix<double>& state_variables, 
        Matrix<double>& forcing);

    // __global__ void AddForcingTerms_kernel(
    //     double* rate_constants, 
    //     int rate_reactants_size, 
    //     double* state_variables, 
    //     double* forcing, 
    //     int matrix_rows, 
    //     int rate_constants_columns, 
    //     int state_forcing_columns,
    //     size_t* number_of_reactants_, 
    //     size_t* accumulated_n_reactants, 
    //     size_t* reactant_ids_,
    //     size_t* number_of_products_, 
    //     size_t* accumulated_n_products, 
    //     size_t* product_ids_, size_t* yields_
    // );

        } // namespace cuda
} // namespace micm
