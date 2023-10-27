// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <chrono>
#include <micm/solver/cuda_lu_decomposition.cuh>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/cuda_param.hpp>
#include <stdexcept>

namespace micm
{
  class CudaLuDecomposition : public LuDecomposition
  {
   public:
    CudaLuDecomposition(){};

    template<typename T, template<class> typename SparseMatrixPolicy>
    requires VectorizableSparse<SparseMatrixPolicy<T>> std::chrono::nanoseconds
    Decompose(const SparseMatrixPolicy<T>& A, SparseMatrixPolicy<T>& L, SparseMatrixPolicy<T>& U)
    const;
  };

  template<typename T, template<class> class SparseMatrixPolicy>
  requires(VectorizableSparse<SparseMatrixPolicy<T>>) std::chrono::nanoseconds
      CudaLuDecomposition::Decompose(const SparseMatrixPolicy<T>& A, SparseMatrixPolicy<T>& L, SparseMatrixPolicy<T>& U)
  const
  {
    CudaSparseMatrixParam sparseMatrix;
    sparseMatrix.A_ = A.AsVector().data();
    sparseMatrix.A_size_ = A.AsVector().size();
    sparseMatrix.L_ = L.AsVector().data();
    sparseMatrix.L_size_ = L.AsVector().size();
    sparseMatrix.U_ = U.AsVector().data();
    sparseMatrix.U_size_ = U.AsVector().size();
    sparseMatrix.n_grids_ = A.size();

    CudaSolverParam solver;
    solver.do_aik_ = do_aik_.data();
    solver.do_aik_size_ = do_aik_.size();
    solver.aik_ = aik_.data();
    solver.aik_size_ = aik_.size();
    solver.do_aki_ = do_aki_.data();
    solver.do_aki_size_ = do_aki_.size();
    solver.aki_ = aki_.data();
    solver.aki_size_ = aki_.size();
    solver.uii_ = uii_.data();
    solver.uii_size_ = uii_.size();

    solver.niLU_ = niLU_.data();
    solver.niLU_size_ = niLU_.size();
    solver.uik_nkj_ = uik_nkj_.data();
    solver.uik_nkj_size_ = uik_nkj_.size();
    solver.lij_ujk_ = lij_ujk_.data();
    solver.lij_ujk_size_ = lij_ujk_.size();
    solver.lki_nkj_ = lki_nkj_.data();
    solver.lki_nkj_size_ = lki_nkj_.size();
    solver.lkj_uji_ = lkj_uji_.data();
    solver.lkj_uji_size_ = lkj_uji_.size();

    // calling kernelSetup function
    return micm::cuda::DecomposeKernelDriver(sparseMatrix, solver);
  }
}  // namespace micm