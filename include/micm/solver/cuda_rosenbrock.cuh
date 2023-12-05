
#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm{
    namespace cuda{
 __global__ void AlphaMinusJacobianKernel(size_t n_grids,
                                          double* d_jacobian,
                                          double* d_jacobian_diagonal_elements,
                                          size_t jacobian_diagonal_elements_size,
                                          const double alpha);

    }//end cuda
}//end micm