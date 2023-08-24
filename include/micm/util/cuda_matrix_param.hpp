#pragma once
#include <micm/util/vector_matrix.hpp>
namespace micm{
    class CUDAMatrixParam{
        public: 
        const double* rate_constants_; 
        template<size_t L>
        CUDAMatrixParam(const VectorMatrix<double, L>& rateConstants){
            rate_constants_ = rateConstants.AsVector().data(); 
        }
    }; 
}