// Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  template<std::size_t L>
  inline JitLuDecompositionDoolittle<L>::JitLuDecompositionDoolittle(JitLuDecompositionDoolittle &&other)
      : LuDecompositionDoolittle(std::move(other)),
        decompose_function_resource_tracker_(std::move(other.decompose_function_resource_tracker_)),
        decompose_function_(std::move(other.decompose_function_))
  {
    other.decompose_function_ = NULL;
  }

  template<std::size_t L>
  inline JitLuDecompositionDoolittle<L> &JitLuDecompositionDoolittle<L>::operator=(JitLuDecompositionDoolittle &&other)
  {
    LuDecompositionDoolittle::operator=(std::move(other));
    decompose_function_resource_tracker_ = std::move(other.decompose_function_resource_tracker_);
    decompose_function_ = std::move(other.decompose_function_);
    other.decompose_function_ = NULL;
    return *this;
  }

  template<std::size_t L>
  inline JitLuDecompositionDoolittle<L>::JitLuDecompositionDoolittle(
      const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix)
  {
    using SparseMatrixPolicy = SparseMatrix<double, SparseMatrixVectorOrdering<L>>;
    decompose_function_ = NULL;
    if (matrix.NumberOfBlocks() > L)
    {
      std::string msg = "JIT functions require the number of grid cells solved together (" +
                        std::to_string(matrix.NumberOfBlocks()) +
                        ") to match the vector dimension template parameter, "
                        "currently: " +
                        std::to_string(L);
      throw std::system_error(make_error_code(MicmJitErrc::InvalidMatrix), msg);
    }
    Initialize<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(
        matrix, typename SparseMatrixPolicy::value_type());
    GenerateDecomposeFunction();
  }

  template<std::size_t L>
  JitLuDecompositionDoolittle<L>::~JitLuDecompositionDoolittle()
  {
    if (decompose_function_ != NULL)
    {
      llvm::ExitOnError exit_on_error;
      exit_on_error(decompose_function_resource_tracker_->remove());
    }
  }

  template<std::size_t L>
  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline JitLuDecompositionDoolittle<L> JitLuDecompositionDoolittle<L>::Create(const SparseMatrixPolicy &matrix)
  {
    static_assert(
        std::is_same_v<SparseMatrixPolicy, LMatrixPolicy>,
        "SparseMatrixPolicy must be the same as LMatrixPolicy for JIT LU decomposition");
    static_assert(
        std::is_same_v<SparseMatrixPolicy, UMatrixPolicy>,
        "SparseMatrixPolicy must be the same as UMatrixPolicy for JIT LU decomposition");
    JitLuDecompositionDoolittle<L> lu_decomp(matrix);
    return lu_decomp;
  }

  template<std::size_t L>
  void JitLuDecompositionDoolittle<L>::GenerateDecomposeFunction()
  {
    std::string function_name = "lu_decompose_" + GenerateRandomString();
    JitFunction func = JitFunction::Create()
                           .SetName(function_name)
                           .SetArguments({ { "A matrix", JitType::DoublePtr },
                                           { "lower matrix", JitType::DoublePtr },
                                           { "upper matrix", JitType::DoublePtr } })
                           .SetReturnType(JitType::Void);
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
        else
        {
          auto loop = func.StartLoop("Uik_eq_zero_loop", 0, L);
          llvm::Value *zero_val = llvm::ConstantFP::get(*(func.context_), llvm::APFloat(0.0));
          llvm::Value *iUf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, uik_nkj->first));
          llvm::Value *U_ptr_index[1];
          U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUf);
          func.SetArrayElement(func.arguments_[2], U_ptr_index, JitType::Double, zero_val);
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
        else
        {
          auto loop = func.StartLoop("Lki_eq_zero_loop", 0, L);
          llvm::Value *zero_val = llvm::ConstantFP::get(*(func.context_), llvm::APFloat(0.0));
          llvm::Value *iLf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, lki_nkj->first));
          llvm::Value *L_ptr_index[1];
          L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLf);
          func.SetArrayElement(func.arguments_[1], L_ptr_index, JitType::Double, zero_val);
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
    decompose_function_ = (FuncPtr)(intptr_t)target.second;
    decompose_function_resource_tracker_ = target.first;
  }

  template<std::size_t L>
  template<class SparseMatrixPolicy>
  void JitLuDecompositionDoolittle<L>::Decompose(const SparseMatrixPolicy &A, auto &lower, auto &upper) const
  {
    decompose_function_(A.AsVector().data(), lower.AsVector().data(), upper.AsVector().data());
    for (size_t block = 0; block < A.NumberOfBlocks(); ++block)
    {
      auto diagonals = upper.DiagonalIndices(block);
    }
  }
}  // namespace micm