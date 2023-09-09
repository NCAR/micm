#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>

namespace micm
{
  namespace cuda
  {
    void DecomposeKernelDriver(
            CUDAMatrixParam& sparseMatrix, 
            CUDASolverParam& solver,
            std::vector<std::pair<std::size_t, std::size_t>>& niLU_,
            std::vector<std::pair<std::size_t, std::size_t>>& uik_nkj_,
            std::vector<std::pair<std::size_t, std::size_t>>& lij_ujk_,
            std::vector<std::pair<std::size_t, std::size_t>>& lki_nkj_,
            std::vector<std::pair<std::size_t, std::size_t>>& lkj_uji_);
  }  // namespace cuda
}  // namespace micm
