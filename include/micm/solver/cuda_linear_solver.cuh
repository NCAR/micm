#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm
{
  namespace cuda
  {
    std::chrono::nanoseconds SolveKernelDriver(
     CudaLinearSolverParam& linearSolver,
     CudaSparseMatrixParam& sparseMatrix, 
     CudaMatrixParam& denseMatrix);
  }  // namespace cuda
}  // namespace micm
