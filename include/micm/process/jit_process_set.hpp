// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process_set.hpp>

namespace micm
{

  /// @brief JIT-compiled solver function calculators for a collection of processes
  ///        The template parameter is the number of grid cells to solve simultaneously
  template<std::size_t L>
  class JitProcessSet : public ProcessSet<VectorMatrix>
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP resource_tracker_;
    void (*forcing_function_)(double*, double*, double*);

   public:
    /// @brief Create a JITed process set calculator for a given set of processes
    /// @param compiler JIT compiler
    /// @param processes Processes to create calculator for
    /// @param state Solver state
    JitProcessSet(
        std::shared_ptr<JitCompiler> compiler,
        const std::vector<Process>& processes,
        const State<VectorMatrix>& state);

    /// @brief Add forcing terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param forcing Forcing terms for each state variable (grid cell, state variable)
    void AddForcingTerms(
        const VectorMatrix<double>& rate_constants,
        const VectorMatrix<double>& state_variables,
        VectorMatrix<double>& forcing) const;
  };

  inline JitProcessSet::JitProcessSet(
      std::shared_ptr<JitCompiler> compiler,
      const std::vector<Process>& processes,
      const State<VectorMatrix>& state)
      : ProcessSet<VectorMatrix>(processes, state),
        compiler_(compiler)
  {
    JitFunction func = JitFunction::create(compiler.get())
                           .name("add_forcing_terms")
                           .arguments({ { "rate constants", JitType::Double },
                                        { "state variables", JitType::Double },
                                        { "forcing", JitType::Double } }),
                .return_type(JitType::Void);
    for (std::size_t i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
    {
      llvm::Value* rc_start = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, i_rxn * L));
      llvm::Value* rc_end = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, i_rxn * L + L));
      llvm::ArrayType* rate_arr = llvm::ArrayType::get(func.GetType(JitType::Int32), ) auto loop =
          func.StartLoop("rate constant loop", rc_start, rc_end);
      llvm::Value* rate = func.GetArrayElement(func.arguments_[0], index_list, micm::JitType::Double);
    }
  }
}  // namespace micm