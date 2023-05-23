#pragma once

#include <vector>

#include <micm/process/process.hpp>
#include <micm/system/system.hpp>

namespace micm
{

  struct State
  {
    System system_;
    std::vector<Process> processes_;
    double temperature_;
    double pressure_;
    std::vector<double> concentrations_;

    /// @brief
    /// @param systems
    State(System system);

    /// @brief
    /// @param system
    /// @param processes
    State(System system, std::vector<Process> processes);

    //
    // TODO: jiwon 5/22 - in progress
    // what is the user's argument type? string or Species?
    //
    /// @brief Set concentrations
    /// @param species_to_concentration
    void set_concentrations(const std::unordered_map<std::string, double>& species_to_concentration);

    /// @brief Update the photolysis rates contained within the processes vector
    void update_photo_rates(double photo_rates[]);

    // std::function<> get_jacobian();
    // std::function<> get_forcing();
  };

  inline State::State(System system)
      : system_(system),
        processes_(),
        temperature_(0),
        pressure_(0),
        concentrations_()
  //
  // TODO: jiwon 5/22
  //
  // concentrations_(nullptr),
  // concentrations_size_(0)
  {
  }

  inline State::State(System system, std::vector<Process> processes)
      : system_(system),
        processes_(processes),
        temperature_(0),
        pressure_(0),
        concentrations_()
  //
  // TODO: jiwon 5/22
  //
  // concentrations_(nullptr),
  // concentrations_size_(0)
  {
  }

  inline void State::set_concentrations(const std::unordered_map<std::string, double>& species_to_concentration)
  {
    concentrations_.reserve(system_.gas_phase_.species_.size());
    for (auto&& species : system_.gas_phase_.species_)
    {
      auto species_ptr = species_to_concentration.find(species.name_);
      if (species_ptr == species_to_concentration.end())
      {
        throw std::runtime_error("Invalid species: " + species.name_);
      }
      concentrations_.push_back(species_ptr->second);
    }
  }

  inline void State::update_photo_rates(double photo_rates[])
  {
    // TODO: do it
  }

}  // namespace micm