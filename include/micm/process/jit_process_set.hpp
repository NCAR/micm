// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/process/process_set.hpp>
#include <micm/util/vector_matrix.hpp>

namespace micm
{

  /// @brief JIT-compiled solver function calculators for a collection of processes
  ///        The template parameter is the number of grid cells to solve simultaneously
  template<std::size_t L>
  class JitProcessSet : public ProcessSet
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP forcing_function_resource_tracker_;
    void (*forcing_function_)(const double *, const double *, double *);

   public:
    /// @brief Create a JITed process set calculator for a given set of processes
    /// @param compiler JIT compiler
    /// @param processes Processes to create calculator for
    /// @param state Solver state
    template<template<class> class MatrixPolicy>
    JitProcessSet(
        std::shared_ptr<JitCompiler> compiler,
        const std::vector<Process> &processes,
        const State<MatrixPolicy> &state);

    ~JitProcessSet();

    /// @brief Add forcing terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param forcing Forcing terms for each state variable (grid cell, state variable)
    template<template<class> class MatrixPolicy>
    void AddForcingTerms(
        const MatrixPolicy<double> &rate_constants,
        const MatrixPolicy<double> &state_variables,
        MatrixPolicy<double> &forcing) const;

   private:
    /// @brief Generate a function to calculate forcing terms
    void GenerateForcingFunction();
  };

  template<std::size_t L>
  template<template<class> class MatrixPolicy>
  inline JitProcessSet<L>::JitProcessSet(
      std::shared_ptr<JitCompiler> compiler,
      const std::vector<Process> &processes,
      const State<MatrixPolicy> &state)
      : ProcessSet(processes, state),
        compiler_(compiler)
  {
    MatrixPolicy<double> test_matrix;
    if (test_matrix.GroupVectorSize() != L)
    {
      std::cerr << "Vector matrix group size invalid for JitProcessSet";
      std::exit(micm::ExitCodes::InvalidMatrixDimension);
    }
    this->GenerateForcingFunction();
  }

  template<std::size_t L>
  void JitProcessSet<L>::GenerateForcingFunction()
  {
    JitFunction func = JitFunction::create(compiler_)
                           .name("add_forcing_terms")
                           .arguments({ { "rate constants", JitType::DoublePtr },
                                        { "state variables", JitType::DoublePtr },
                                        { "forcing", JitType::DoublePtr } })
                           .return_type(JitType::Void);
    llvm::Type *double_type = func.GetType(JitType::Double);
    llvm::Value *zero = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 0));
    llvm::Type *rate_array_type = llvm::ArrayType::get(double_type, L);
    llvm::AllocaInst *rate_array = func.builder_->CreateAlloca(
        rate_array_type, llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 1)), "rate_array");

    auto react_ids = reactant_ids_.begin();
    auto prod_ids = product_ids_.begin();
    auto yields = yields_.begin();
    for (std::size_t i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
    {
      llvm::Value *rc_start = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, i_rxn * L));

      // save rate constant in rate array for each grid cell
      auto loop = func.StartLoop("rate constant", 0, L);
      llvm::Value *ptr_index[1];
      ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, rc_start);
      llvm::Value *rate_const = func.GetArrayElement(func.arguments_[0], ptr_index, JitType::Double);
      llvm::Value *array_index[2];
      array_index[0] = zero;
      array_index[1] = loop.index_;
      llvm::Value *rate_ptr = func.builder_->CreateInBoundsGEP(rate_array_type, rate_array, array_index);
      func.builder_->CreateStore(rate_const, rate_ptr);
      func.EndLoop(loop);

      // rates[i_cell] *= reactant_concentration for each reactant
      for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
      {
        loop = func.StartLoop("rate calc", 0, L);
        llvm::Value *react_id = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, react_ids[i_react] * L));
        ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, react_id);
        llvm::Value *react_conc = func.GetArrayElement(func.arguments_[1], ptr_index, JitType::Double);
        array_index[1] = loop.index_;
        rate_ptr = func.builder_->CreateInBoundsGEP(rate_array_type, rate_array, array_index);
        llvm::Value *rate = func.builder_->CreateLoad(double_type, rate_ptr, "rate");
        rate = func.builder_->CreateFMul(rate, react_conc, "rate");
        func.builder_->CreateStore(rate, rate_ptr);
        func.EndLoop(loop);
      }

      // set forcing for each reactant f[i_react][i_cell] -= rate[i_cell]
      for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
      {
        loop = func.StartLoop("reactant forcing", 0, L);
        llvm::Value *react_id = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, react_ids[i_react] * L));
        ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, react_id);
        llvm::Value *react_forcing_ptr = func.builder_->CreateGEP(double_type, func.arguments_[2].ptr_, ptr_index);
        llvm::Value *react_forcing = func.builder_->CreateLoad(double_type, react_forcing_ptr, "forcing");
        array_index[1] = loop.index_;
        rate_ptr = func.builder_->CreateInBoundsGEP(rate_array_type, rate_array, array_index);
        llvm::Value *rate = func.builder_->CreateLoad(double_type, rate_ptr, "rate");
        react_forcing = func.builder_->CreateFSub(react_forcing, rate, "forcing");
        func.builder_->CreateStore(react_forcing, react_forcing_ptr);
        func.EndLoop(loop);
      }

      // set forcing for each product f[i_prod][i_cell] += yield * rate[i_cell]
      for (std::size_t i_prod = 0; i_prod < number_of_products_[i_rxn]; ++i_prod)
      {
        loop = func.StartLoop("product forcing", 0, L);
        llvm::Value *prod_id = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, prod_ids[i_prod] * L));
        ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, prod_id);
        llvm::Value *prod_forcing_ptr = func.builder_->CreateGEP(double_type, func.arguments_[2].ptr_, ptr_index);
        llvm::Value *prod_forcing = func.builder_->CreateLoad(double_type, prod_forcing_ptr, "forcing");
        array_index[1] = loop.index_;
        rate_ptr = func.builder_->CreateInBoundsGEP(rate_array_type, rate_array, array_index);
        llvm::Value *rate = func.builder_->CreateLoad(double_type, rate_ptr, "rate");
        llvm::Value *yield = llvm::ConstantFP::get(*(func.context_), llvm::APFloat(yields[i_prod]));
        rate = func.builder_->CreateFMul(rate, yield, "rate_yield");
        prod_forcing = func.builder_->CreateFAdd(prod_forcing, rate, "forcing");
        func.builder_->CreateStore(prod_forcing, prod_forcing_ptr);
        func.EndLoop(loop);
      }
      react_ids += number_of_reactants_[i_rxn];
      prod_ids += number_of_products_[i_rxn];
      yields += number_of_products_[i_rxn];
    }
    func.builder_->CreateRetVoid();

    auto target = func.Generate();
    forcing_function_ = (void (*)(const double *, const double *, double *))(intptr_t)target.second;
    forcing_function_resource_tracker_ = target.first;
  }

  template<std::size_t L>
  JitProcessSet<L>::~JitProcessSet()
  {
    if (forcing_function_resource_tracker_)
    {
      llvm::ExitOnError exit_on_error;
      exit_on_error(forcing_function_resource_tracker_->remove());
    }
  }

  template<std::size_t L>
  template<template<class> class MatrixPolicy>
  void JitProcessSet<L>::AddForcingTerms(
      const MatrixPolicy<double> &rate_constants,
      const MatrixPolicy<double> &state_variables,
      MatrixPolicy<double> &forcing) const
  {
    forcing_function_(rate_constants.AsVector().data(), state_variables.AsVector().data(), forcing.AsVector().data());
  }

}  // namespace micm