
#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm{
    namespace cuda{
  void AlphaMinusJacobianDriver(
                        CudaSparseMatrixParam& sparseMatrix,
                        std::vector<double> jacobian_diagonal_elements, 
                        const double alpha);

    }//end cuda
}//end micm