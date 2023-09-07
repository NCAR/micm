// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix optimized with JIT-compilation
  ///
  /// See LuDecomposition class description for algorithm details
  template<std::size_t L>
  class JitLuDecomposition : public LuDecomposition
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP decompose_function_resource_tracker_;
    void (*decompose_function_)(const double *, double *, double *);

   public:
    /// @brief Create a JITed LU decomposer for a given sparse matrix structure
    /// @param compiler JIT compiler
    /// @param matrix Sparse matrix to create LU decomposer for
    JitLuDecomposition(
        std::shared_ptr<JitCompiler> compiler,
        const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix);

    ~JitLuDecomposition();

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @param lower The lower triangular matrix created by decomposition
    /// @param upper The upper triangular matrix created by decomposition
    template<typename T, template<class> class SparseMatrixPolicy>
    void Decompose(const SparseMatrixPolicy<T> &A, SparseMatrixPolicy<T> &lower, SparseMatrixPolicy<T> &upper) const;

   private:
    /// @brief Generates a function to perform the LU decomposition for a specific matrix sparsity structure
    void GenerateDecomposeFunction();
  };

  template<std::size_t L>
  inline JitLuDecomposition<L>::JitLuDecomposition(
      std::shared_ptr<JitCompiler> compiler,
      const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix)
      : LuDecomposition(matrix),
        compiler_(compiler)
  {
    decompose_function_ = NULL;
    assert(matrix.size() <= L && "Jit LU Decomposition matrix size mismatch");
    GenerateDecomposeFunction();
  }

  template<std::size_t L>
  JitLuDecomposition<L>::~JitLuDecomposition()
  {
    if (decompose_function_ != NULL)
    {
      llvm::ExitOnError exit_on_error;
      exit_on_error(decompose_function_resource_tracker_->remove());
    }
  }

  template<std::size_t L>
  void JitLuDecomposition<L>::GenerateDecomposeFunction()
  {
    JitFunction func = JitFunction::create(compiler_)
                           .name("lu_decompose")
                           .arguments({ { "A matrix", JitType::DoublePtr },
                                        { "lower matrix", JitType::DoublePtr },
                                        { "upper matrix", JitType::DoublePtr } })
                           .return_type(JitType::Void);
    llvm::Type *double_type = func.GetType(JitType::Double);
    llvm::Value *zero = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 0));
    auto do_aik = do_aik_.begin();
    auto aik = aik_.begin();
    auto uik_nkj = uik_nkj_.begin();
    auto lij_ujk = lij_ujk_.begin();
    auto do_aki = do_aki_.begin();
    auto aki = aki_.begin();
    auto lki_nkj = lki_nkj_.begin();
    auto lkj_uji = lkj_uji_.begin();
    auto uii = uii_.begin();
    for (auto &inLU : niLU_)
    {
      // Upper triangular matrix
      for (std::size_t iU = 0; iU < inLU.second; ++iU)
      {
        if (*(do_aik++))
        {
          auto loop = func.StartLoop("Uik_eq_Aik_loop", 0, L);
          llvm::Value *iAf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, *(aik++)));
          llvm::Value *A_ptr_index[1];
          A_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iAf);
          llvm::Value *A_val = func.GetArrayElement(func.arguments_[0], A_ptr_index, JitType::Double);
          llvm::Value *iUf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, uik_nkj->first));
          llvm::Value *U_ptr_index[1];
          U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUf);
          func.SetArrayElement(func.arguments_[2], U_ptr_index, JitType::Double, A_val);
          func.EndLoop(loop);
        }
        for (std::size_t ikj = 0; ikj < uik_nkj->second; ++ikj)
        {
          auto loop = func.StartLoop("Uik_seq_Lij_Ujk_loop", 0, L);
          llvm::Value *iLf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lij_ujk->first));
          llvm::Value *L_ptr_index[1];
          L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLf);
          llvm::Value *L_val = func.GetArrayElement(func.arguments_[1], L_ptr_index, JitType::Double);
          llvm::Value *iUf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lij_ujk->second));
          llvm::Value *U_ptr_index[1];
          U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUf);
          llvm::Value *U_val = func.GetArrayElement(func.arguments_[2], U_ptr_index, JitType::Double);
          llvm::Value *prod = func.builder_->CreateFMul(L_val, U_val, "Lij_mul_Ujk");
          iUf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, uik_nkj->first));
          U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUf);
          llvm::Value *U_ptr = func.builder_->CreateGEP(double_type, func.arguments_[2].ptr_, U_ptr_index);
          U_val = func.builder_->CreateLoad(double_type, U_ptr);
          U_val = func.builder_->CreateFSub(U_val, prod, "Uik_seq_Lij_Ujk");
          func.builder_->CreateStore(U_val, U_ptr);
          func.EndLoop(loop);
          ++lij_ujk;
        }
        ++uik_nkj;
      }
      // Lower triangular matrix
      {
        auto loop = func.StartLoop("Lki_eq_1_loop", 0, L);
        llvm::Value *L_val = llvm::ConstantFP::get(*(func.context_), llvm::APFloat(1.0));
        llvm::Value *iLf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, (lki_nkj++)->first));
        llvm::Value *L_ptr_index[1];
        L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLf);
        func.SetArrayElement(func.arguments_[1], L_ptr_index, JitType::Double, L_val);
        func.EndLoop(loop);
      }
      for (std::size_t iL = 0; iL < inLU.first; ++iL)
      {
        if (*(do_aki++))
        {
          auto loop = func.StartLoop("Lki_eq_Aki", 0, L);
          llvm::Value *iAf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, *(aki++)));
          llvm::Value *A_ptr_index[1];
          A_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iAf);
          llvm::Value *A_val = func.GetArrayElement(func.arguments_[0], A_ptr_index, JitType::Double);
          llvm::Value *iLf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lki_nkj->first));
          llvm::Value *L_ptr_index[1];
          L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLf);
          func.SetArrayElement(func.arguments_[1], L_ptr_index, JitType::Double, A_val);
          func.EndLoop(loop);
        }
        for (std::size_t ikj = 0; ikj < lki_nkj->second; ++ikj)
        {
          auto loop = func.StartLoop("Lki_seq_Lkj_Uji_loop", 0, L);
          llvm::Value *iLf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lkj_uji->first));
          llvm::Value *L_ptr_index[1];
          L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLf);
          llvm::Value *L_val = func.GetArrayElement(func.arguments_[1], L_ptr_index, JitType::Double);
          llvm::Value *iUf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lkj_uji->second));
          llvm::Value *U_ptr_index[1];
          U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUf);
          llvm::Value *U_val = func.GetArrayElement(func.arguments_[2], U_ptr_index, JitType::Double);
          llvm::Value *prod = func.builder_->CreateFMul(L_val, U_val, "Lkj_mul_Uji");
          iLf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lki_nkj->first));
          L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLf);
          llvm::Value *L_ptr = func.builder_->CreateGEP(double_type, func.arguments_[1].ptr_, L_ptr_index);
          L_val = func.builder_->CreateLoad(double_type, L_ptr);
          L_val = func.builder_->CreateFSub(L_val, prod, "Lki_seq_Lkj_Uji");
          func.builder_->CreateStore(L_val, L_ptr);
          func.EndLoop(loop);
          ++lkj_uji;
        }
        {
          auto loop = func.StartLoop("Lki_deq_Uii_loop", 0, L);
          llvm::Value *iUf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, *(uii++)));
          llvm::Value *U_ptr_index[1];
          U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUf);
          llvm::Value *U_val = func.GetArrayElement(func.arguments_[2], U_ptr_index, JitType::Double);
          llvm::Value *iLf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lki_nkj->first));
          llvm::Value *L_ptr_index[1];
          L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLf);
          llvm::Value *L_ptr = func.builder_->CreateGEP(double_type, func.arguments_[1].ptr_, L_ptr_index);
          llvm::Value *L_val = func.builder_->CreateLoad(double_type, L_ptr);
          L_val = func.builder_->CreateFDiv(L_val, U_val, "Lki_deq_Uii");
          func.builder_->CreateStore(L_val, L_ptr);
          func.EndLoop(loop);
          ++lki_nkj;
        }
      }
    }

    func.builder_->CreateRetVoid();

    auto target = func.Generate();
    decompose_function_ = (void (*)(const double *, double *, double *))(intptr_t)target.second;
    decompose_function_resource_tracker_ = target.first;
  }

  template<std::size_t L>
  template<typename T, template<class> class SparseMatrixPolicy>
  void JitLuDecomposition<L>::Decompose(
      const SparseMatrixPolicy<T> &A,
      SparseMatrixPolicy<T> &lower,
      SparseMatrixPolicy<T> &upper) const
  {
    decompose_function_(A.AsVector().data(), lower.AsVector().data(), upper.AsVector().data());
  }
}  // namespace micm