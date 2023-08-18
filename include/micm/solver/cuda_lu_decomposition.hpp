#pragma once
#include<micm/solver/lu_decomposition.hpp>
#include<thrust/device_vector.h> 
#ifdef USE_CUDA
#include <micm/solver/cuda_de_composition.cuh>
#endif 

namespace micm{
    class CUDALuDecomposition: public LuDecomposition{
        public: 
        CUDALuDecomposition(); 
        /// @brief Construct an LU decomposition algorithm for a given sparse matrix
        /// @param matrix Sparse matrix
        template<typename T, typename OrderingPolicy>
        CUDALuDecomposition(const SparseMatrix<T, OrderingPolicy>& matrix); 
    
    #ifdef USE_CUDA
        template<typename T, template<class> typename SparseMatrixPolicy>
        requires VectorizableSparse<SparseMatrixPolicy<T>>
        void Decompose(
            const SparseMatrixPolicy<T>&A, 
            SparseMatrixPolicy<T>& L, 
            SparseMatrixPolicy<T>& U) const; 
    #endif
    
    inline CUDALuDecomposition::CUDALuDecomposition(){};

    template<typename T, typename OrderingPolicy>
    inline CUDALuDecomposition::CUDALuDecomposition(const SparseMatrix<T, OrderingPolicy>& matrix){
        std::size_t n = matrix[0].size();
        auto LU = GetLUMatrices(matrix, T{});
        const auto& L_row_start = LU.first.RowStartVector();
        const auto& L_row_ids = LU.first.RowIdsVector();
        const auto& U_row_start = LU.second.RowStartVector();
        const auto& U_row_ids = LU.second.RowIdsVector();
        for (std::size_t i = 0; i < matrix[0].size(); ++i)
        {
        std::pair<std::size_t, std::size_t> iLU(0, 0);
        // Upper triangular matrix
        for (std::size_t k = i; k < n; ++k)
        {
            std::size_t nkj = 0;
            for (std::size_t j_id = L_row_start[i]; j_id < L_row_start[i + 1]; ++j_id)
            {
            std::size_t j = L_row_ids[j_id];
            if (j >= i)
                break;
            if (LU.second.IsZero(j, k))
                continue;
            ++nkj;
            lij_ujk_.push_back(std::make_pair(LU.first.VectorIndex(0, i, j), LU.second.VectorIndex(0, j, k)));
            }
            if (matrix.IsZero(i, k))
            {
            if (nkj == 0 && k != i)
                continue;
            do_aik_.push_back(false);
            }
            else
            {
            do_aik_.push_back(true);
            aik_.push_back(matrix.VectorIndex(0, i, k));
            }
            uik_nkj_.push_back(std::make_pair(LU.second.VectorIndex(0, i, k), nkj));
            ++(iLU.second);
        }
        // Lower triangular matrix
        lki_nkj_.push_back(std::make_pair(LU.first.VectorIndex(0, i, i), 0));
        for (std::size_t k = i + 1; k < n; ++k)
        {
            std::size_t nkj = 0;
            for (std::size_t j_id = L_row_start[k]; j_id < L_row_start[k + 1]; ++j_id)
            {
            std::size_t j = L_row_ids[j_id];
            if (j >= i)
                break;
            if (LU.second.IsZero(j, i))
                continue;
            ++nkj;
            lkj_uji_.push_back(std::make_pair(LU.first.VectorIndex(0, k, j), LU.second.VectorIndex(0, j, i)));
            }
            if (matrix.IsZero(k, i))
            {
            if (nkj == 0)
                continue;
            do_aki_.push_back(false);
            }
            else
            {
            do_aki_.push_back(true);
            aki_.push_back(matrix.VectorIndex(0, k, i));
            }
            uii_.push_back(LU.second.VectorIndex(0, i, i));
            lki_nkj_.push_back(std::make_pair(LU.first.VectorIndex(0, k, i), nkj));
            ++(iLU.first);
        }
        niLU_.push_back(iLU);
    }
  } 

    template<typename T, template<class> class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy<T>>) 
    void CUDALuDecomposition::Decompose(
        const SparseMatrixPolicy<T>& A, 
        SparseMatrixPolicy<T>& L, 
        SparseMatrixPolicy<T>& U) const
    {
        thrust::device_vector d_niLU<std::pair<size_t,size_t>> = niLU_; 
        thrust::device_vector d_uik_nkj<std::pair<size_t,size_t>> = uik_nkj_;
        thrust::device_vector d_lij_ujk<std::pair<size_t, size_t>> = lij_ujk_; 
        thrust::device_vector d_lki_nkj<std::pair<size_t, size_t>> = lki_nkj_; 
        thrust::device_vector d_lkj_uji<std::pair<size_t, size_t>> = lkj_uji_;
        
        //calling kernelSetup function
        
    }


    } //end class CUDALuDecomposition
}//end micm 