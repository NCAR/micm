#pragma once
#include <vector> 

  struct CUDAProcessSetParam{
    const size_t* number_of_reactants; 
    const size_t* reactant_ids; 
    size_t reactant_ids_size; 
    const size_t* number_of_products; 
    const size_t* product_ids; 
    size_t product_ids_size; 
    double* yields; 
    size_t yields_size; 
    const size_t* jacobian_flat_ids; 
    size_t jacobian_flat_ids_size; 
}; 
namespace micm{
    class CUDAMatrixParam{
        public: 
        const double* rate_constants_; 
        const double* state_variables_; 
        double* forcing_; 
        double* jacobian_; 
        size_t n_grids_; 
        size_t n_reactions_; 
        size_t n_species_; 
        size_t jacobian_size_; 
        CUDAMatrixParam(){};
    }; //end class
}//end micm
