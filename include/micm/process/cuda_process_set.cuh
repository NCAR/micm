#pragma once
#include <micm/util/cuda_matrix_param.hpp>
namespace micm
{
  namespace cuda
  {
    std::chrono::nanoseconds AddForcingTermsKernelDriver(
        micm::CUDAMatrixParam& matrixParam,
        const size_t* number_of_reactants,
        const size_t* reactant_ids,
        size_t reactant_ids_size,
        const size_t* number_of_products,
        const size_t* product_ids,
        size_t product_ids_size,
        const double* yields,
        size_t yields_size); 

    std::chrono::nanoseconds AddJacobianTermsKernelDriver(
        const double* rate_constants,
        const double* state_variables,
        size_t n_grids,
        size_t n_reactions,
        size_t n_species,
        double* jacobian,
        size_t jacobian_size,
        const size_t* number_of_reactants, 
        const size_t* reactant_ids,
        size_t reactant_ids_size,
        const size_t* number_of_products,
        const double* yields,
        size_t yields_size,
        const size_t* jacobian_flat_ids,
        size_t jacobian_flat_ids_size);
  }  // namespace cuda
}  // namespace micm
