// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/chemical_reaction.hpp>
#include <micm/process/phase_transfer_process.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/yield.hpp>
#include <micm/util/error.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  class Process
  {
   public:
    using ProcessVariant = std::variant<ChemicalReaction, PhaseTransferProcess>;

    ProcessVariant process_;

    template<typename T>
      requires std::same_as<std::decay_t<T>, ChemicalReaction> || std::same_as<std::decay_t<T>, PhaseTransferProcess>
    Process(T&& process)
        : process_(std::forward<T>(process))
    {
    }

    /// @brief TODO - Temporary wrapper for rate constant calculation
    ///        Calls ChemicalReaction::CalculateRateConstants for all ChemicalReaction processes
    ///        issue - https://github.com/NCAR/micm/issues/812
    template<
        class DenseMatrixPolicy,
        class SparseMatrixPolicy,
        class LuDecompositionPolicy,
        class LMatrixPolicy,
        class UMatrixPolicy>
    static void CalculateRateConstants(
        const std::vector<Process>& processes,
        State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state)
    {
      // Collect ChemicalReaction objects
      std::vector<ChemicalReaction> reactions;
      for (const auto& process : processes)
      {
        if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
        {
          reactions.push_back(*reaction);
        }
        // PhaseTransferProcess support can be added here in the future
      }

      ChemicalReaction::CalculateRateConstants(reactions, state);
    }
  };

}  // namespace micm
