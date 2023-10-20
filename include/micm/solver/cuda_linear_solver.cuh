#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm
{
  namespace cuda
  {
    void SolveKernelDriver(
     CudaLinearSolverParam& linearSolver,
     CudaSparseMatrixParam& sparseMatrix, 
     CudaMatrixParam& denseMatrix);
  }  // namespace cuda
}  // namespace micm
