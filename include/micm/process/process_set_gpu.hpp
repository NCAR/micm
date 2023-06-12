#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <vector>

namespace micm {

    void AddForcingTerms_kernelSetup(
        const Matrix<double>& rate_constants,
        const Matrix<double>& state_variables,
        Matrix<double>& forcing,
        std::vector<std::size_t> number_of_reactants_,
        std::vector<std::size_t> reactant_ids_,
        std::vector<std::size_t> number_of_products_,
        std::vector<std::size_t> product_ids_,
        std::vector<std::size_t> yields_
    );

    __global__ void AddForcingTerms_kernel(
        std::vector<double>* rate_constants,
        std::vector<double> state_variables,
        std::vector<double> forcing,
        int matrix_rows,
        int rate_constants_columns,
        int state_forcing_columns,
        std::vector<std::size_t> number_of_reactants_,
        std::vector<std::size_t> accumulated_n_reactants,
        std::vector<std::size_t> reactant_ids_,
        std::vector<std::size_t> number_of_products_,
        std::vector<std::size_t> accumulated_n_products,
        std::vector<std::size_t> product_ids_,
        std::vector<std::size_t> yields_
    );

} // namespace micm
