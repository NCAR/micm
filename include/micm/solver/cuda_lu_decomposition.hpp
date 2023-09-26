#pragma once
#include<micm/solver/lu_decomposition.hpp>
#include<micm/util/cuda_param.hpp>
#include <stdexcept>
#ifdef USE_CUDA
#include <micm/solver/cuda_lu_decomposition.cuh>
#endif 

#ifdef USE_CUDA
namespace micm{
    class CUDALuDecomposition: public LuDecomposition{
    public: 
    CUDALuDecomposition(){}; 
    
    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<typename T, typename OrderingPolicy>
    CUDALuDecomposition(const SparseMatrix<T, OrderingPolicy>& matrix); 
    
    template<typename T, template<class> typename SparseMatrixPolicy>
    requires VectorizableSparse<SparseMatrixPolicy<T>>
    void Decompose(
        const SparseMatrixPolicy<T>&A, 
        SparseMatrixPolicy<T>& L, 
        SparseMatrixPolicy<T>& U) const; 
    }; 

    template<typename T, template<class> class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy<T>>) 
    void CUDALuDecomposition::Decompose(
        const SparseMatrixPolicy<T>& A, 
        SparseMatrixPolicy<T>& L, 
        SparseMatrixPolicy<T>& U) const
    {
        CUDASparseMatrixParam sparseMatrix; 
        sparseMatrix.A = A.AsVector().data(); 
        sparseMatrix.A_size = A.AsVector().size(); 
        sparseMatrix.L = L.AsVector().data(); 
        sparseMatrix.L_size = L.AsVector().size(); 
        sparseMatrix.U = U.AsVector().data(); 
        sparseMatrix.U_size = U.AsVector().size(); 
        sparseMatrix.n_grids = A.size(); 
        
        CUDASolverParam solver;    
        solver.do_aik = do_aik_.data(); 
        solver.do_aik_size = do_aik_.size(); 
        solver.aik = aik_.data(); 
        solver.aik_size = aik_.size(); 
        solver.do_aki = do_aki_.data(); 
        solver.do_aki_size = do_aki_.size(); 
        solver.aki = aki_.data(); 
        solver.aki_size = aki_.size(); 
        solver.uii = uii_.data(); 
        solver.uii_size = uii_.size(); 
        
        solver.niLU = niLU_.data(); 
        solver.niLU_size = niLU_.size(); 
        solver.uik_nkj = uik_nkj_.data(); 
        solver.uik_nkj_size = uik_nkj_.size(); 
        solver.lij_ujk = lij_ujk_.data(); 
        solver.lij_ujk_size = lij_ujk_.size(); 
        solver.lki_nkj = lki_nkj_.data(); 
        solver.lki_nkj_size = lki_nkj_.size(); 
        solver.lkj_uji =  lkj_uji_.data(); 
        solver.lkj_uji_size = lkj_uji_.size(); 

        //calling kernelSetup function
        micm::cuda::DecomposeKernelDriver(sparseMatrix, solver); 
    }
}//end micm
#endif