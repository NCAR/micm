// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/chemical_reaction.hpp>
#include <micm/process/phase_transfer_process.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/error.hpp>

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
  };

}  // namespace micm
