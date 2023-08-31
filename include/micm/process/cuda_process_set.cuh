#pragma once
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    std::chrono::nanoseconds AddForcingTermsKernelDriver(
        CUDAMatrixParam& matrixParam,
        CUDAProcessSetParam& processSet); 

    std::chrono::nanoseconds AddJacobianTermsKernelDriver(
        CUDAMatrixParam& matrixParam, 
        CUDAProcessSetParam& processSet);
  }  // namespace cuda
}  // namespace micm
