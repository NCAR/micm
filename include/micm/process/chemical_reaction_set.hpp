// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/chemical_reaction.hpp>
#include <micm/process/process.hpp>
#include <micm/system/species.hpp>
#include <micm/system/yield.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{
  CreateProcessInfoFrom(
    const ChemicalReaction* process,
    const std::pair<std::string, std::size_t>& independent_variable,
    std::size_t i_process,
    const std::map<std::string, std::size_t>& variable_map,
    std::vector<ProcessInfo>& process_infos)
  {
    MICM_PROFILE_FUNCTION();

    // Does this process depend on this species as a reactant?
    bool matches_independent = false;
    for (const auto& ind_reactant : process->reactants_)
    {
      // If the current species appears in the process reactants, 
      // it's an independent variable for this process.
      if (ind_reactant.name_ == independent_variable.first)
      {
        matches_independent = true;
        break;
      }
    }

    if (!matches_independent)
      return;

    ProcessInfo info;
    info.process_id_ = i_process;
    info.independent_id_ = independent_variable.second;

    // Collect contributors: all reactants (except the first occurrence of the independent variable), and all products
    bool found_independent = false;
    for (const auto& reactant : process->reactants_)
    {
      if (reactant.IsParameterized())
        continue;  // Skip parameterized reactants

      auto idx_it = variable_map.find(reactant.name_);
      if (idx_it == variable_map.end())
        throw std::system_error(make_error_code(MicmProcessSetErrc::ReactantDoesNotExist), reactant.name_);

      // Skip the first occurrence of the independent variable
      if (idx_it->second == independent_variable.second && !found_independent)
      {
        found_independent = true;
        continue;
      }
      info.contributors_.push_back({ JacobianContributionType::Reactant, idx_it->second });
    }

    for (const auto& product : process.products_)
    {
      if (product.first.IsParameterized())
        continue;  // Skip parameterized products

      auto idx_it = variable_map.find(product.first.name_);
      if (idx_it == variable_map.end())
        throw std::system_error(make_error_code(MicmProcessSetErrc::ProductDoesNotExist), product.first.name_);

      info.contributors_.push_back({ JacobianContributionType::Product, idx_it->second, product.second });
    }

    process_infos.push_back(info);
  }
  //   }
  // };


}  // namespace micm