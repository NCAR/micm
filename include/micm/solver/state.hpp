#pragma once

#include <cstddef>
#include <vector>

namespace micm
{

  struct State
  {
    double temperature_;
    double pressure_;
    double air_density_;
    std::vector<double> concentrations_;
    std::vector<double> custom_rate_parameters_;
    std::vector<double> rate_constants_;

    /// @brief
    State();

    /// @brief
    /// @param state_size The number of System state variables
    /// @param custom_parameters_size The number of custom rate parameters
    /// @param process_size The number of processes to store rate constants for
    State(const std::size_t state_size, const std::size_t custom_parameters_size, const std::size_t process_size);
  };

  inline State::State()
      : temperature_(0),
        pressure_(0),
        air_density_(1),
        concentrations_(),
        custom_rate_parameters_(),
        rate_constants_()
  {
  }

  inline State::State(const std::size_t state_size, const std::size_t custom_parameters_size, const std::size_t process_size)
      : temperature_(0),
        pressure_(0),
        air_density_(1),
        concentrations_(state_size, 0.0),
        custom_rate_parameters_(custom_parameters_size, 0.0),
        rate_constants_(process_size, 0.0)
  {
  }
}  // namespace micm