  #pragma once

  #include <vector>

  #include <micm/system/system.hpp>
  #include <micm/process/process.hpp>
  
  namespace micm {

    struct State  {
      System system_;
      std::vector<Process> processes_;
      double temperature_;
      double pressure_; 
      std::vector<double> concentrations_;

      /// @brief 
      State();

      /// @brief 
      /// @param system 
      /// @param processes 
      State(System system, std::vector<Process> processes);

      /// @brief Update the photolysis rates contained within the processes vector
      void update_photo_rates();

      // std::function<> get_jacobian();
      // std::function<> get_forcing();
    };

    inline State::State()
      : system_()
      , processes_()
      , temperature_(0)
      , pressure_(0)
      , concentrations_()
    {
    }

    inline State::State(System system, std::vector<Process> processes)
      : system_(system)
      , processes_(processes)
      , temperature_(0)
      , pressure_(0)
      , concentrations_()
    {
    }

    inline void State::update_photo_rates()
    {
      // TODO: do it
    }
  }