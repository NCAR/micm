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
      std::vector<std::pair<std::size_t, std::size_t>>& niLU;
      std::vector<std::pair<std::size_t, std::size_t>>& uik_nkj;
      std::vector<std::pair<std::size_t, std::size_t>>& lij_ujk;
      std::vector<std::pair<std::size_t, std::size_t>>& lki_nkj;
      std::vector<std::pair<std::size_t, std::size_t>>& lkj_uji ; 
      char* do_aik;
      size_t do_aik_size; 
      const size_t* aik;
      size_t aik_size; 
      char* do_aki;
      size_t do_aki_size; 
      const size_t* aki; 
      size_t aki_size; 
      const size_t* uii; 
      size_t uii_size; 
  };
 //different matrix data grouped in struct passing to kernel driver function 
    struct CUDAMatrixParam{
     const double* rate_constants; 
     const double* state_variables; 
     double* forcing; 
     size_t n_grids; 
     size_t n_reactions; 
     size_t n_species; 
}; 
//sparseMatrix data grouped in struct passing to kernel driver function 
struct CUDASparseMatrixParam{
   double* jacobian; 
   size_t jacobian_size; 
   const double* A; 
   size_t A_size; 
   double* L; 
   size_t L_size; 
   double* U;
   size_t U_size; 
};
#endif