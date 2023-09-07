#pragma once
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
   void DecomposeKernelDriver(
            CUDAMatrixParam& sparseMatrix, 
            CUDASolverParam& solver);
  }  // namespace cuda
}  // namespace micm
