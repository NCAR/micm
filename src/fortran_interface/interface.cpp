#include "interface.hpp"

#include <iostream>

namespace micm
{

  void solver(
      double state[],
      int64_t state_size,
      int64_t
          timestep)  // NOLINT(misc-unused-parameters,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
  {
    std::cout << "state size: " << state_size << std::endl;
    std::cout << "timestep: " << timestep << std::endl;

    for (int64_t i{}; i < state_size; ++i)
    {
      std::cout << "state[" << i << "]: " << state[i] << std::endl;
    }
  }

  FuncPtr get_solver(
      char filepath[])  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
  {
    std::cout << "file path: " << filepath << "\n";
    return &solver;
  }

}  // namespace micm