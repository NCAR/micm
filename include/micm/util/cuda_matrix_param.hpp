#pragma once
#include <micm/util/vector_matrix.hpp>
namespace micm{
    class CUDAMatrixParam{
        public: 
        double* rate_constants; 
        template<size_t L>
        CUDAMatrixParam(const VectorMatrix<double, L>& rateConstants){
            rate_constants = rateConstants.AsVector().data(); 
        }
    }; 
}