#include<thrust/device_vector.h> 
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

  struct CUDASolverParam{
      thrust::device_vector d_niLU<thrust::pair<size_t,size_t>>;
      thrust::device_vector d_uik_nkj<thrust::pair<size_t,size_t>>; 
      thrust::device_vector d_lij_ujk<thrust::pair<size_t, size_t>>;
      thrust::device_vector d_lki_nkj<thrust::pair<size_t, size_t>>; 
      thrust::device_vector d_lkj_uji<thrust::pair<size_t, size_t>>;
      const bool* do_aik;
      size_t do_aik_size; 
      const size_t* aik;
      size_t aik_size; 
      const bool* do_aki;
      size_t do_aki_size; 
      const size_t* aki; 
      size_t aki_size; 
      const size_t* uii; 
      size_t uii_size; 
  }
 //different matrix data grouped in struct passing to kernel driver function 
    struct CUDAMatrixParam{
     const double* rate_constants_; 
     const double* state_variables_; 
     double* forcing_; 
     double* jacobian_; 
     const double* A; 
     double* L; 
     double* U; 

     size_t n_grids_; 
     size_t n_reactions_; 
     size_t n_species_; 
     size_t jacobian_size_; 
     size_t A_size; 
     size_t L_size;
     size_t U_size; 
}; 
#endif