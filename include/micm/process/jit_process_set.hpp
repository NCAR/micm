// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/process/process_set.hpp>
#include <micm/util/random_string.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

namespace micm
{

  /// @brief JIT-compiled solver function calculators for a collection of processes
  ///        The template parameter is the number of grid cells to solve simultaneously
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class JitProcessSet : public ProcessSet
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP forcing_function_resource_tracker_;
    void (*forcing_function_)(const double *, const double *, double *);
    llvm::orc::ResourceTrackerSP jacobian_function_resource_tracker_;
    void (*jacobian_function_)(const double *, const double *, double *);

   public:
    JitProcessSet(const JitProcessSet &) = delete;
    JitProcessSet &operator=(const JitProcessSet &) = delete;
    JitProcessSet(JitProcessSet &&);
    JitProcessSet &operator=(JitProcessSet &&);

    JitProcessSet() = default;

    /// @brief Create a JITed process set calculator for a given set of processes
    /// @param compiler JIT compiler
    /// @param processes Processes to create calculator for
    /// @param variable_map A mapping of species names to concentration index
    JitProcessSet(
        std::shared_ptr<JitCompiler> compiler,
        const std::vector<Process> &processes,
        const std::map<std::string, std::size_t> &variable_map);

    ~JitProcessSet();

    /// @brief Sets the indices for each non-zero Jacobian element in the underlying vector
    /// @param matrix The sparse matrix used for the Jacobian
    template<typename OrderingPolicy>
    void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy> &matrix);

    /// @brief Adds forcing terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param forcing Forcing terms for each state variable (grid cell, state variable)
    template<template<class> class MatrixPolicy>
    void AddForcingTerms(
        const MatrixPolicy<double> &rate_constants,
        const MatrixPolicy<double> &state_variables,
        MatrixPolicy<double> &forcing) const;

    /// @brief Subtracts Jacobian terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param jacobian Jacobian matrix for the system (grid cell, dependent variable, independent variable)
    template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
    void SubtractJacobianTerms(
        const MatrixPolicy<double> &rate_constants,
        const MatrixPolicy<double> &state_variables,
        SparseMatrixPolicy<double> &jacobian) const;

   private:
    /// @brief Generates a function to calculate forcing terms
    /// @param matrix The matrix that will hold the forcing terms
    void GenerateForcingFunction();
    /// @brief Generate a function to calculate Jacobian terms
    /// @param matrix The sparse matrix that will hold the Jacobian
    void GenerateJacobianFunction(const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix);
  };

  template<std::size_t L>
  inline JitProcessSet<L>::JitProcessSet(JitProcessSet &&other)
      : ProcessSet(std::move(other)),
        compiler_(std::move(other.compiler_)),
        forcing_function_resource_tracker_(std::move(other.forcing_function_resource_tracker_)),
        forcing_function_(std::move(other.forcing_function_)),
        jacobian_function_resource_tracker_(std::move(other.jacobian_function_resource_tracker_)),
        jacobian_function_(std::move(other.jacobian_function_))
  {
    other.forcing_function_ = NULL;
    other.jacobian_function_ = NULL;
  }

  template<std::size_t L>
  inline JitProcessSet<L> &JitProcessSet<L>::operator=(JitProcessSet &&other)
  {
    ProcessSet::operator=(std::move(other));
    compiler_ = std::move(other.compiler_);
    forcing_function_resource_tracker_ = std::move(other.forcing_function_resource_tracker_);
    forcing_function_ = std::move(other.forcing_function_);
    jacobian_function_resource_tracker_ = std::move(other.jacobian_function_resource_tracker_);
    jacobian_function_ = std::move(other.jacobian_function_);
    other.forcing_function_ = NULL;
    other.jacobian_function_ = NULL;
    return *this;
  }

  template<std::size_t L>
  inline JitProcessSet<L>::JitProcessSet(
      std::shared_ptr<JitCompiler> compiler,
      const std::vector<Process> &processes,
      const std::map<std::string, std::size_t> &variable_map)
      : ProcessSet(processes, variable_map),
        compiler_(compiler)
  {
    forcing_function_ = NULL;
    jacobian_function_ = NULL;
    this->GenerateForcingFunction();
  }

  template<std::size_t L>
  void JitProcessSet<L>::GenerateForcingFunction()
  {
    std::string function_name = "add_forcing_terms_" + generate_random_string();
    JitFunction func = JitFunction::create(compiler_)
                           .name(function_name)
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
  template<typename OrderingPolicy>
  inline void JitProcessSet<L>::SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy> &matrix)
  {
    ProcessSet::SetJacobianFlatIds(matrix);
    GenerateJacobianFunction(matrix);
  }

  template<std::size_t L>
  void JitProcessSet<L>::GenerateJacobianFunction(const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix)
  {
    std::string function_name = "subtract_jacobian_terms_" + generate_random_string();
    JitFunction func = JitFunction::create(compiler_)
                           .name(function_name)
                           .arguments({ { "rate constants", JitType::DoublePtr },
                                        { "state variables", JitType::DoublePtr },
                                        { "jacobian", JitType::DoublePtr } })
                           .return_type(JitType::Void);
    llvm::Type *double_type = func.GetType(JitType::Double);
    llvm::Value *zero = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 0));
    llvm::Type *rate_array_type = llvm::ArrayType::get(double_type, L);
    llvm::AllocaInst *rate_array = func.builder_->CreateAlloca(
        rate_array_type, llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 1)), "d_rate_d_ind_array");

    auto react_ids = reactant_ids_.begin();
    auto yields = yields_.begin();
    auto flat_id = jacobian_flat_ids_.begin();
    for (std::size_t i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
    {
      llvm::Value *rc_start = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, i_rxn * L));
      for (std::size_t i_ind = 0; i_ind < number_of_reactants_[i_rxn]; ++i_ind)
      {
        // save rate constant in d_rate_d_ind for each grid cell
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

        // d_rate_d_ind[i_cell] *= reactant_concentration for each reactant except ind
        for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
        {
          if (i_react == i_ind)
            continue;
          loop = func.StartLoop("d_rate_d_ind calc", 0, L);
          llvm::Value *react_id = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, react_ids[i_react] * L));
          ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, react_id);
          llvm::Value *react_conc = func.GetArrayElement(func.arguments_[1], ptr_index, JitType::Double);
          array_index[1] = loop.index_;
          rate_ptr = func.builder_->CreateInBoundsGEP(rate_array_type, rate_array, array_index);
          llvm::Value *rate = func.builder_->CreateLoad(double_type, rate_ptr, "d_rate_d_ind");
          rate = func.builder_->CreateFMul(rate, react_conc, "d_rate_d_ind");
          func.builder_->CreateStore(rate, rate_ptr);
          func.EndLoop(loop);
        }

        // set jacobian terms for each reactant jac[i_react][i_ind][i_cell] += d_rate_d_ind[i_cell]
        for (std::size_t i_dep = 0; i_dep < number_of_reactants_[i_rxn]; ++i_dep)
        {
          loop = func.StartLoop("reactant term", 0, L);
          llvm::Value *dep_id = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, *(flat_id++)));
          ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, dep_id);
          llvm::Value *dep_jac_ptr = func.builder_->CreateGEP(double_type, func.arguments_[2].ptr_, ptr_index);
          llvm::Value *dep_jac = func.builder_->CreateLoad(double_type, dep_jac_ptr, "reactant jacobian term");
          array_index[1] = loop.index_;
          rate_ptr = func.builder_->CreateInBoundsGEP(rate_array_type, rate_array, array_index);
          llvm::Value *rate = func.builder_->CreateLoad(double_type, rate_ptr, "d_rate_d_ind");
          dep_jac = func.builder_->CreateFAdd(dep_jac, rate, "reactant jacobian term");
          func.builder_->CreateStore(dep_jac, dep_jac_ptr);
          func.EndLoop(loop);
        }

        // set jacobian terms for each product jac[i_prod][i_ind][i_cell] -= yield * d_rate_d_ind[i_cell]
        for (std::size_t i_dep = 0; i_dep < number_of_products_[i_rxn]; ++i_dep)
        {
          loop = func.StartLoop("product term", 0, L);
          llvm::Value *dep_id = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, *(flat_id++)));
          ptr_index[0] = func.builder_->CreateNSWAdd(loop.index_, dep_id);
          llvm::Value *dep_jac_ptr = func.builder_->CreateGEP(double_type, func.arguments_[2].ptr_, ptr_index);
          llvm::Value *dep_jac = func.builder_->CreateLoad(double_type, dep_jac_ptr, "product jacobian term");
          array_index[1] = loop.index_;
          rate_ptr = func.builder_->CreateInBoundsGEP(rate_array_type, rate_array, array_index);
          llvm::Value *rate = func.builder_->CreateLoad(double_type, rate_ptr, "d_rate_d_ind");
          llvm::Value *yield = llvm::ConstantFP::get(*(func.context_), llvm::APFloat(yields[i_dep]));
          rate = func.builder_->CreateFMul(rate, yield, "product yield");
          dep_jac = func.builder_->CreateFSub(dep_jac, rate, "product jacobian term");
          func.builder_->CreateStore(dep_jac, dep_jac_ptr);
          func.EndLoop(loop);
        }
      }
      react_ids += number_of_reactants_[i_rxn];
      yields += number_of_products_[i_rxn];
    }
    func.builder_->CreateRetVoid();

    auto target = func.Generate();
    jacobian_function_ = (void (*)(const double *, const double *, double *))(intptr_t)target.second;
    jacobian_function_resource_tracker_ = target.first;
  }

  template<std::size_t L>
  JitProcessSet<L>::~JitProcessSet()
  {
    if (forcing_function_ != NULL)
    {
      llvm::ExitOnError exit_on_error;
      exit_on_error(forcing_function_resource_tracker_->remove());
    }
    if (jacobian_function_ != NULL)
    {
      llvm::ExitOnError exit_on_error;
      exit_on_error(jacobian_function_resource_tracker_->remove());
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

  template<std::size_t L>
  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  void JitProcessSet<L>::SubtractJacobianTerms(
      const MatrixPolicy<double> &rate_constants,
      const MatrixPolicy<double> &state_variables,
      SparseMatrixPolicy<double> &jacobian) const
  {
    jacobian_function_(rate_constants.AsVector().data(), state_variables.AsVector().data(), jacobian.AsVector().data());
  }

}  // namespace micm