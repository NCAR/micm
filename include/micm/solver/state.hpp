  #pragma once

  #include <cstddef>
  #include <vector>

  namespace micm {

    struct State  {
      double temperature_;
      double pressure_;
      double air_density_;
      std::vector<double> concentrations_;
      std::vector<double> custom_rate_parameters_;

      /// @brief 
      State();

      /// @brief 
      /// @param state_size The number of System state variables 
      /// @param rate_parameter_size The number of user-defined process parameters
      State(const std::size_t state_size, const std::size_t rate_parameter_size);
    };

    inline State::State()
      : temperature_(0)
      , pressure_(0)
      , air_density_(1)
      , concentrations_()
      , custom_rate_parameters_()
    {
    }

    inline State::State(const std::size_t state_size, const std::size_t rate_parameter_size)
      : temperature_(0)
      , pressure_(0)
      , air_density_(1)
      , concentrations_(state_size, 0.0)
      , custom_rate_parameters_(rate_parameter_size, 0.0)
    {
    }
  }