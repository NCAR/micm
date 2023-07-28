#pragma once

namespace micm
{
  namespace cuda
  {
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
        int yields_size);
<<<<<<< HEAD:include/micm/process/process_set_cuda.cuh

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
        size_t reactants_ids_size, 
        const size_t* number_of_products, 
        const size_t* product_ids,
        size_t product_ids_size, 
        const double* yields,
        size_t yields_size,
        const size_t* jacobian_flat_ids,
        size_t jacobian_flat_ids_size);
        } // namespace cuda
} // namespace micm
=======
  }  // namespace cuda
}  // namespace micm
>>>>>>> main:include/micm/process/cuda_process_set.cuh
