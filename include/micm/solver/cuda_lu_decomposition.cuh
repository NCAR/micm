#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm
{
  namespace cuda
  {
    void DecomposeKernelDriver(
            CUDASparseMatrixParam& sparseMatrix, 
            CUDASolverParam& solver);
  }  // namespace cuda
}  // namespace micm