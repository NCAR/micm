#pragma once
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    std::chrono::nanoseconds AddForcingTermsKernelDriver(
        CudaMatrixParam& matrix,
        CudaProcessSetParam& processSet); 

    std::chrono::nanoseconds AddJacobianTermsKernelDriver(
        CudaMatrixParam& matrix, 
        CudaSparseMatrixParam& sparseMatrix, 
        CudaProcessSetParam& processSet);
  }  // namespace cuda
}  // namespace micm
