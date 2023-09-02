#ifndef CUDA_PARAM_HPP
#define CUDA_PARAM_HPP
  //member data of class CUDAProcessSet grouped in struct passing to kernel driver function 
  struct CUDAProcessSetParam{
    const size_t* number_of_reactants; 
    const size_t* reactant_ids; 
    size_t reactant_ids_size; 
    const size_t* number_of_products; 
    const size_t* product_ids; 
    size_t product_ids_size; 
    const double* yields; 
    size_t yields_size; 
    const size_t* jacobian_flat_ids; 
    size_t jacobian_flat_ids_size; 
}; 
//different matrix data grouped in struct passing to kernel driver function 
    struct CUDAMatrixParam{
     const double* rate_constants; 
     const double* state_variables; 
     double* forcing; 
     double* jacobian; 
     size_t n_grids; 
     size_t n_reactions; 
     size_t n_species; 
     size_t jacobian_size; 
}; 
//sparseMatrix data grouped in struct passing to kernel driver function 
struct CUDASparseMatrixParam{
   double* jacobian; 
   size_t jacobian_size; 
};
#endif