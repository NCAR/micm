// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
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

enum class MicmProcessErrc
{
  TooManyReactantsForSurfaceReaction = MICM_PROCESS_ERROR_CODE_TOO_MANY_REACTANTS_FOR_SURFACE_REACTION
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmProcessErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class MicmProcessErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_PROCESS;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmProcessErrc>(ev))
      {
        case MicmProcessErrc::TooManyReactantsForSurfaceReaction: return "A surface reaction can only have one reactant";
        default: return "Unknown error";
      }
    }
  };

  const MicmProcessErrorCategory MICM_PROCESS_ERROR{};
}  // namespace

inline std::error_code make_error_code(MicmProcessErrc e)
{
  return { static_cast<int>(e), MICM_PROCESS_ERROR };
}

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
