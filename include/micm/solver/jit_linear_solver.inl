// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::JitLinearSolver(JitLinearSolver &&other)
      : LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>(std::move(other)),
        compiler_(std::move(other.compiler_)),
        solve_function_resource_tracker_(std::move(other.solve_function_resource_tracker_)),
        solve_function_(std::move(other.solve_function_))
  {
    other.solve_function_ = NULL;
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy> &
  JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::operator=(JitLinearSolver &&other)
  {
    LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::operator=(std::move(other));
    compiler_ = std::move(other.compiler_);
    solve_function_resource_tracker_ = std::move(other.solve_function_resource_tracker_);
    solve_function_ = std::move(other.solve_function_);
    other.solve_function_ = NULL;
    return *this;
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::JitLinearSolver(
      std::shared_ptr<JitCompiler> compiler,
      const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix,
      double initial_value)
      : LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>(
            matrix,
            initial_value,
            [&](const SparseMatrixPolicy<double> &m) -> LuDecompositionPolicy
            { return LuDecompositionPolicy(compiler, m); }),
        compiler_(compiler)
  {
    solve_function_ = NULL;
    if (matrix.size() != L || matrix.GroupVectorSize() != L)
    {
      throw std::runtime_error("Invalid matrix for JitLinearSolver. Check the the VectorMatrix template parameters.");
    }
    GenerateSolveFunction();
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::~JitLinearSolver()
  {
    if (solve_function_ != NULL)
    {
      llvm::ExitOnError exit_on_error;
      exit_on_error(solve_function_resource_tracker_->remove());
    }
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& lower_matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& upper_matrix,
        bool& is_singular)
  {
    LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::Factor(matrix, lower_matrix, upper_matrix, is_singular);
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& lower_matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& upper_matrix)
  {
    LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::Factor(matrix, lower_matrix, upper_matrix);
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  template<template<class> class MatrixPolicy>
  inline void JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      const MatrixPolicy<double> &b,
      MatrixPolicy<double> &x,
      SparseMatrixPolicy<double> &lower_matrix,
      SparseMatrixPolicy<double> &upper_matrix)
  {
    solve_function_(
        b.AsVector().data(),
        x.AsVector().data(),
        lower_matrix.AsVector().data(),
        upper_matrix.AsVector().data());
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void JitLinearSolver<L, SparseMatrixPolicy, LuDecompositionPolicy>::GenerateSolveFunction()
  {
    std::string function_name = "linear_solve_" + generate_random_string();
    JitFunction func = JitFunction::create(compiler_)
                           .name(function_name)
                           .arguments({ { "b", JitType::DoublePtr },
                                        { "x", JitType::DoublePtr },
                                        { "L", JitType::DoublePtr },
                                        { "U", JitType::DoublePtr } })
                           .return_type(JitType::Void);
    llvm::Type *double_type = func.GetType(JitType::Double);
    auto Lij_yj = LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::Lij_yj_.begin();
    auto Uij_xj = LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::Uij_xj_.begin();
    std::size_t offset = 0;
    for (auto &nLij_Lii : LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::nLij_Lii_)
    {
      // the x vector is used for y values to conserve memory
      {
        auto loop = func.StartLoop("yi_eq_bi_loop", 0, L);
        llvm::Value *ibf = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, offset));
        llvm::Value *ptr_index[1];
        ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, ibf);
        llvm::Value *b_val = func.GetArrayElement(func.arguments_[0], ptr_index, JitType::Double);
        func.SetArrayElement(func.arguments_[1], ptr_index, JitType::Double, b_val);
        func.EndLoop(loop);
      }
      for (std::size_t i = 0; i < nLij_Lii.first; ++i)
      {
        auto loop = func.StartLoop("yi_seq_Lij_yj_loop", 0, L);
        llvm::Value *iyj = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, (*Lij_yj).second * L));
        llvm::Value *y_ptr_index[1];
        y_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iyj);
        llvm::Value *yj_val = func.GetArrayElement(func.arguments_[1], y_ptr_index, JitType::Double);
        llvm::Value *iLij = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, (*Lij_yj).first));
        llvm::Value *L_ptr_index[1];
        L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLij);
        llvm::Value *L_val = func.GetArrayElement(func.arguments_[2], L_ptr_index, JitType::Double);
        yj_val = func.builder_->CreateFMul(L_val, yj_val, "Lij_mul_yj");
        llvm::Value *iyi = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, offset));
        y_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iyi);
        llvm::Value *yi_ptr = func.builder_->CreateGEP(double_type, func.arguments_[1].ptr_, y_ptr_index);
        llvm::Value *yi_val = func.builder_->CreateLoad(double_type, yi_ptr);
        yi_val = func.builder_->CreateFSub(yi_val, yj_val, "yi_seq_Lij_yj");
        func.builder_->CreateStore(yi_val, yi_ptr);
        func.EndLoop(loop);
        ++Lij_yj;
      }
      {
        auto loop = func.StartLoop("yi_deq_Lii_loop", 0, L);
        llvm::Value *iLii = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, nLij_Lii.second));
        llvm::Value *L_ptr_index[1];
        L_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iLii);
        llvm::Value *L_val = func.GetArrayElement(func.arguments_[2], L_ptr_index, JitType::Double);
        llvm::Value *iyi = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, offset));
        llvm::Value *y_ptr_index[1];
        y_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iyi);
        llvm::Value *y_ptr = func.builder_->CreateGEP(double_type, func.arguments_[1].ptr_, y_ptr_index);
        llvm::Value *y_val = func.builder_->CreateLoad(double_type, y_ptr);
        y_val = func.builder_->CreateFDiv(y_val, L_val, "yi_deq_Lii");
        func.builder_->CreateStore(y_val, y_ptr);
        func.EndLoop(loop);
      }
      offset += L;
    }
    offset = L * LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::nUij_Uii_.size();
    for (auto &nUij_Uii : LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>::nUij_Uii_)
    {
      offset -= L;
      for (std::size_t i = 0; i < nUij_Uii.first; ++i)
      {
        auto loop = func.StartLoop("xi_seq_Uij_xj_loop", 0, L);
        llvm::Value *ixj = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, (*Uij_xj).second * L));
        llvm::Value *x_ptr_index[1];
        x_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, ixj);
        llvm::Value *xj_val = func.GetArrayElement(func.arguments_[1], x_ptr_index, JitType::Double);
        llvm::Value *iUij = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, (*Uij_xj).first));
        llvm::Value *U_ptr_index[1];
        U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUij);
        llvm::Value *U_val = func.GetArrayElement(func.arguments_[3], U_ptr_index, JitType::Double);
        xj_val = func.builder_->CreateFMul(xj_val, U_val, "Uij_mul_xj");
        llvm::Value *ixi = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, offset));
        x_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, ixi);
        llvm::Value *xi_ptr = func.builder_->CreateGEP(double_type, func.arguments_[1].ptr_, x_ptr_index);
        llvm::Value *xi_val = func.builder_->CreateLoad(double_type, xi_ptr);
        xi_val = func.builder_->CreateFSub(xi_val, xj_val, "xi_seq_Uij_xj");
        func.builder_->CreateStore(xi_val, xi_ptr);
        func.EndLoop(loop);
        ++Uij_xj;
      }
      {
        auto loop = func.StartLoop("xi_deq_Uii_loop", 0, L);
        llvm::Value *iUii = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, nUij_Uii.second));
        llvm::Value *U_ptr_index[1];
        U_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, iUii);
        llvm::Value *U_val = func.GetArrayElement(func.arguments_[3], U_ptr_index, JitType::Double);
        llvm::Value *ixi = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, offset));
        llvm::Value *x_ptr_index[1];
        x_ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, ixi);
        llvm::Value *x_ptr = func.builder_->CreateGEP(double_type, func.arguments_[1].ptr_, x_ptr_index);
        llvm::Value *x_val = func.builder_->CreateLoad(double_type, x_ptr);
        x_val = func.builder_->CreateFDiv(x_val, U_val, "xi_deq_Uii");
        func.builder_->CreateStore(x_val, x_ptr);
        func.EndLoop(loop);
      }
    }

    func.builder_->CreateRetVoid();

    auto target = func.Generate();
    solve_function_ = (void (*)(const double *, double *, const double *, const double *))(intptr_t)target.second;
    solve_function_resource_tracker_ = target.first;
  }

}  // namespace micm