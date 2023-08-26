#pragma once
#include <micm/util/cuda_matrix_param.hpp>
#include <micm/process/cuda_process_set.hpp>
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
        micm::CUDAMatrixParam& matrixParam, 
        CUDAProcessSetParam& processSetParam, 
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
