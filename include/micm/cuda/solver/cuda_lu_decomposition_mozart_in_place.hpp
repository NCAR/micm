// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>

namespace micm
{
  /// This CudaLuDecompositionMozartInPlace class inherits everything from the base class "LuDecompositionMozartInPlace"
  class CudaLuDecompositionMozartInPlace : public LuDecompositionMozartInPlace
  {
   public:
    /// This is an instance of struct "LuDecomposeMozartInPlaceParam" that holds
    ///   the constant data of "CudaLuDecompositionMozartInPlace" class on the device
    LuDecomposeMozartInPlaceParam devstruct_;

    /// This is the default constructor, taking no arguments;
    CudaLuDecompositionMozartInPlace(){};

    CudaLuDecompositionMozartInPlace(const CudaLuDecompositionMozartInPlace&) = delete;
    CudaLuDecompositionMozartInPlace& operator=(const CudaLuDecompositionMozartInPlace&) = delete;
    CudaLuDecompositionMozartInPlace(CudaLuDecompositionMozartInPlace&& other)
        : LuDecompositionMozartInPlace(std::move(other))
    {
      std::swap(this->devstruct_, other.devstruct_);
    };

    CudaLuDecompositionMozartInPlace& operator=(CudaLuDecompositionMozartInPlace&& other)
    {
      LuDecompositionMozartInPlace::operator=(std::move(other));
      std::swap(this->devstruct_, other.devstruct_);
      return *this;
    };

    /// This is the overloaded constructor that takes one argument called "matrix";
    /// We need to specify the type (e.g., double, int, etc) and
    ///   ordering (e.g., vector-stored, non-vector-stored, etc) of the "matrix";
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    CudaLuDecompositionMozartInPlace(const SparseMatrixPolicy& matrix)
    {
      Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());

      /// Convert host-side tuple/pair vectors to uint32_t arrays with packed pairs
      const std::size_t n = this->aii_nji_nki_.size();

      // Extract SoA from aii_nji_nki_ tuples
      std::vector<uint32_t> h_aii(n), h_nji(n), h_nki(n);
      for (std::size_t i = 0; i < n; ++i)
      {
        h_aii[i] = static_cast<uint32_t>(std::get<0>(this->aii_nji_nki_[i]));
        h_nji[i] = static_cast<uint32_t>(std::get<1>(this->aii_nji_nki_[i]));
        h_nki[i] = static_cast<uint32_t>(std::get<2>(this->aii_nji_nki_[i]));
      }

      // Convert aji_ from size_t to uint32_t
      std::vector<uint32_t> h_aji(this->aji_.size());
      for (std::size_t i = 0; i < this->aji_.size(); ++i)
        h_aji[i] = static_cast<uint32_t>(this->aji_[i]);

      // Pack (aik, njk) pairs interleaved: [aik0, njk0, aik1, njk1, ...]
      // Kernel loads as uint2 for a single 8-byte transaction per pair
      std::vector<uint32_t> h_aik_njk_packed(this->aik_njk_.size() * 2);
      for (std::size_t i = 0; i < this->aik_njk_.size(); ++i)
      {
        h_aik_njk_packed[i * 2] = static_cast<uint32_t>(this->aik_njk_[i].first);
        h_aik_njk_packed[i * 2 + 1] = static_cast<uint32_t>(this->aik_njk_[i].second);
      }

      // Pack (ajk, aji_update) pairs interleaved: [ajk0, aji0, ajk1, aji1, ...]
      // Kernel loads as uint2 for a single 8-byte transaction per pair
      std::vector<uint32_t> h_ajk_aji_packed(this->ajk_aji_.size() * 2);
      for (std::size_t i = 0; i < this->ajk_aji_.size(); ++i)
      {
        h_ajk_aji_packed[i * 2] = static_cast<uint32_t>(this->ajk_aji_[i].first);
        h_ajk_aji_packed[i * 2 + 1] = static_cast<uint32_t>(this->ajk_aji_[i].second);
      }

      /// Populate the host struct
      LuDecomposeMozartInPlaceParam hoststruct;
      hoststruct.aii_ = h_aii.data();
      hoststruct.nji_ = h_nji.data();
      hoststruct.nki_ = h_nki.data();
      hoststruct.aji_ = h_aji.data();
      hoststruct.aik_njk_packed_ = h_aik_njk_packed.data();
      hoststruct.ajk_aji_packed_ = h_ajk_aji_packed.data();
      hoststruct.n_ = static_cast<uint32_t>(n);
      hoststruct.aji_size_ = static_cast<uint32_t>(this->aji_.size());
      hoststruct.aik_njk_size_ = static_cast<uint32_t>(this->aik_njk_.size());
      hoststruct.ajk_aji_size_ = static_cast<uint32_t>(this->ajk_aji_.size());

      /// Create the ALU matrix with all the fill-ins for the non-zero values
      auto ALU = GetLUMatrix<SparseMatrixPolicy>(matrix, 0, true);
      hoststruct.number_of_non_zeros_ = static_cast<uint32_t>(ALU.GroupSize() / SparseMatrixPolicy::GroupVectorSize());

      // Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    };

    /// This is destructor that will free the device memory of
    ///   the constant data from the class "CudaLuDecompositionMozartInPlace"
    ~CudaLuDecompositionMozartInPlace()
    {
      /// Free the device memory allocated by the members of "devstruct_"
      micm::cuda::FreeConstData(this->devstruct_);
    };

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static CudaLuDecompositionMozartInPlace Create(const SparseMatrixPolicy& matrix)
    {
      CudaLuDecompositionMozartInPlace lu_decomp(matrix);
      return lu_decomp;
    }

    /// @brief This is the function to perform an LU decomposition on a given A matrix on the GPU
    /// @param ALU Sparse matrix to decompose (will be overwritten with L and U matrices)
    template<class SparseMatrixPolicy>
      requires(CudaMatrix<SparseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
    void Decompose(SparseMatrixPolicy& ALU) const;
  };

  template<class SparseMatrixPolicy>
    requires(CudaMatrix<SparseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
  void CudaLuDecompositionMozartInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    auto ALU_param = ALU.AsDeviceParam();
    micm::cuda::DecomposeKernelDriver(ALU_param, this->devstruct_);
  }
}  // end of namespace micm
