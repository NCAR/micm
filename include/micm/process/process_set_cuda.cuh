#pragma once


#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <vector>

namespace micm {
    namespace cuda {

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
        Matrix<double>& forcing);



        } // namespace cuda
} // namespace micm
