#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm{
    namespace cuda{
  void AlphaMinusJacobianDriver(
                        CudaSparseMatrixParam& sparseMatrix,
                      
                        double alpha);

    }//end cuda
}//end micm