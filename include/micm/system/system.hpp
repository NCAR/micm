/* Copyright (C) 2023-2024 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <functional>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace micm
{
  struct SystemParameters
  {
    Phase gas_phase_{};
    std::unordered_map<std::string, Phase> phases_{};
  };

  /**
   * @brief A `System` holds all physical information that represents a grid cell.
   *
   */
  class System
  {
   public:
    /// @brief The gas phase is a micm::Phase and determines what species are present.
    Phase gas_phase_;
    /// @brief This is a catchall for anything that is not the gas phase.
    std::unordered_map<std::string, Phase> phases_;

   public:
    /// @brief Default constructor
    System()
        : gas_phase_(),
          phases_()
    {
    }

    System(const Phase& gas_phase, const std::unordered_map<std::string, Phase>& phases)
        : gas_phase_(gas_phase),
          phases_(phases)
    {
    }

    System(Phase&& gas_phase, std::unordered_map<std::string, Phase>&& phases)
        : gas_phase_(std::move(gas_phase)),
          phases_(std::move(phases))
    {
    }

    System(const System& other)
        : gas_phase_(other.gas_phase_),
          phases_(other.phases_)
    {
    }

    System(System&& other)
        : gas_phase_(std::move(other.gas_phase_)),
          phases_(std::move(other.phases_))
    {
    }

    System(const SystemParameters& parameters)
        : gas_phase_(parameters.gas_phase_),
          phases_(parameters.phases_)
    {
    }

    System& operator=(const System& other)
    {
      gas_phase_ = other.gas_phase_;
      phases_ = other.phases_;
      return *this;
    }

    /// @brief Returns the number of doubles required to store the system state
    size_t StateSize() const;

    /// @brief Returns a set of unique species names
    /// @return vector of unique state variable names
    std::vector<std::string> UniqueNames() const;

    /// @brief Returns a set of unique species names
    /// @param f Function used to apply specific order to unique names
    /// @return vector of unique state variable names
    std::vector<std::string> UniqueNames(
        const std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> f) const;
  };

  inline size_t System::StateSize() const
  {
    size_t state_size = gas_phase_.StateSize();
    for (const auto& phase : phases_)
    {
      state_size += phase.second.StateSize();
    }
    return state_size;
  }

  inline std::vector<std::string> System::UniqueNames() const
  {
    return UniqueNames(nullptr);
  }

  inline std::vector<std::string> System::UniqueNames(
      const std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> f) const
  {
    std::vector<std::string> names = gas_phase_.UniqueNames();
    for (auto& phase : phases_)
    {
      for (auto& species_name : phase.second.UniqueNames())
        names.push_back(phase.first + "." + species_name);
    }
    if (f)
    {
      const auto orig_names = names;
      for (std::size_t i = 0; i < orig_names.size(); ++i)
        names[i] = f(orig_names, i);
    }
    return names;
  }

}  // namespace micm
