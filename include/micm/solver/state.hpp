  #pragma once

  #include <vector>
  
  namespace micm {

    struct State  {
      double temperature;
      double pressure; 
      std::vector<double> concentrations;

      // std::function<> get_jacobian();
      // std::function<> get_forcing();
    };

  }