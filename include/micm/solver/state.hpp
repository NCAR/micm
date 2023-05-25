#pragma once

#include <cstddef>
#include <vector>

namespace micm
{

  struct StateParameters
  {
    std::size_t number_of_grid_cells_{ 1 };
    std::size_t number_of_state_variables_;
    std::size_t number_of_custom_parameters_{ 0 };
    std::size_t number_of_rate_constants_{ 0 };
  };

  struct Conditions{
    double temperature_{0.0};
    double pressure_{0.0};
    double air_density_{1.0};
  };

  struct State
  {
    std::vector<Conditions> conditions_;
    std::vector<std::vector<double>> concentrations_;
    std::vector<std::vector<double>> custom_rate_parameters_;
    std::vector<std::vector<double>> rate_constants_;

    /// @brief
    State();

    /// @brief
    /// @param state_size The number of System state variables
    /// @param custom_parameters_size The number of custom rate parameters
    /// @param process_size The number of processes to store rate constants for
    State(const std::size_t state_size, const std::size_t custom_parameters_size, const std::size_t process_size);

    /// @brief
    /// @param parameters State dimension information
    State(const StateParameters parameters);
  };

  inline State::State()
      : conditions_(),
        concentrations_(),
        custom_rate_parameters_(),
        rate_constants_()
  {
  }

  inline State::State(const std::size_t state_size, const std::size_t custom_parameters_size, const std::size_t process_size)
      : conditions_(1),
        concentrations_(1, std::vector<double>( state_size, 0.0 )),
        custom_rate_parameters_(1, std::vector<double>( custom_parameters_size, 0.0 )),
        rate_constants_(1, std::vector<double>( process_size, 0.0 ))
  {
  }

  inline State::State(const StateParameters parameters)
      : conditions_(parameters.number_of_grid_cells_),
        concentrations_(parameters.number_of_grid_cells_, std::vector<double>( parameters.number_of_state_variables_, 0.0 )),
        custom_rate_parameters_(parameters.number_of_grid_cells_, std::vector<double>( parameters.number_of_custom_parameters_, 0.0 )),
        rate_constants_(parameters.number_of_grid_cells_, std::vector<double>( parameters.number_of_rate_constants_, 0.0 ))
  {
  }
}  // namespace micm