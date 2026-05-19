// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/chemical_reaction.hpp>

#include <utility>

namespace micm
{

  class Process
  {
   public:
    ChemicalReaction process_;

    Process(const ChemicalReaction& process) 
        : process_(process)
    {
    }

    Process(ChemicalReaction&& process)
         : process_(std::move(process))
    {
    }
  };

}  // namespace micm
