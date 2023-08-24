#pragma once
#include <vector> 
namespace micm{
    class CUDAMatrixParam{
        public: 
        const double* rate_constants_; 
        CUDAMatrixParam(const vector<double> rate_constants){
            rate_constants_ = rate_constants.AsVector().data(); 
        }
    }; 
}