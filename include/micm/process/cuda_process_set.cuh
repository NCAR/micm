#pragma once
#include <micm/util/cuda_matrix_param.hpp>

namespace micm
{
  namespace cuda
  {
    std::chrono::nanoseconds AddForcingTermsKernelDriver(
        CUDAMatrixParam& matrixParam,
        const size_t* number_of_reactants,
        const size_t* reactant_ids,
        size_t reactant_ids_size,
        const size_t* number_of_products,
        const size_t* product_ids,
        size_t product_ids_size,
        const double* yields,
        size_t yields_size); 

    std::chrono::nanoseconds AddJacobianTermsKernelDriver(
        CUDAMatrixParam& matrixParam, 
        CUDAProcessSetParam& processSet);
  }  // namespace cuda
}  // namespace micm
