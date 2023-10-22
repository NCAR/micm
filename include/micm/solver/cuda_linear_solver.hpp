#pragma once
#include <micm/solver/linear_solver.hpp> 
#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/util/cuda_param.hpp>

#ifdef USE_CUDA
#include <micm/solver/cuda_linear_solver.cuh>
namespace micm{
    template<typename T, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy = CudaLuDecomposition>
    class CudaLinearSolver: public LinearSolver<T, SparseMatrixPolicy, LuDecomposition> {
    public:
        //constructor
        CudaLinearSolver(){};

        CudaLinearSolver(const SparseMatrixPolicy<T>& matrix, T initial_value): LinearSolver<T, SparseMatrixPolicy, LuDecomposition> (matrix, initial_value){};
      
        template<template<class> class MatrixPolicy> 
        requires(VectorizableDense<MatrixPolicy<T>> || VectorizableSparse<SparseMatrixPolicy<T>>)
        void Solve(const MatrixPolicy<T>&b, MatrixPolicy<T>& x)
        {
            CudaLinearSolverParam linearSolver;
            CudaSparseMatrixParam sparseMatrix; 
            CudaMatrixParam denseMatrix;
            linearSolver.nLij_Lii_ = nLij_Lii_.data(); 
            linearSolver.nLij_Lii_size_ = nLij_Lii_.size(); 
            linearSolver.Lij_yj_= Lij_yj_.data();
            linearSolver.Lij_yj_size_ = Lij_yj_.size();
            linearSolver.nUij_Uii_ = nUij_Uii_.data();
            linearSolver.nUij_Uii_size_ = nUij_Uii_.size();
            linearSolver.Uij_xj_ = Uij_xj_.data(); 
            linearSolver.Uij_xj_size_ = Uij_xj_.size(); 
            sparseMatrix.lower_matrix_ = lower_matrix_.AsVector().data(); 
            sparseMatrix.lower_matrix_size_ = lower_matrix_.AsVector().size();
            sparseMatrix.upper_matrix_ = upper_matrix_.AsVector().data(); 
            sparseMatrix.upper_matrix_size_ = upper_matrix_.AsVector().size(); 
            denseMatrix.b_ = b.AsVector().data(); 
            denseMatrix.x_ = x.AsVector().data();
            denseMatrix.n_grids_ = b.size(); //number of grids
            denseMatrix.b_column_counts_ = b[0].size(); 
            denseMatrix.x_column_counts_= x[0].size(); 
            
            //calling kernel driver
            micm::cuda::SolveKernelDriver(linearSolver, sparseMatrix, denseMatrix);
        };
    };
}//end micm
#endif