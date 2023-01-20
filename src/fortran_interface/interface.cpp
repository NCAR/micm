#include "interface.hpp"

#include <iostream>

namespace micm
{

  void solver(double state[], uint64_t state_size, uint64_t timestep)
  {  // NOLINT(misc-unused-parameters,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
    std::cout << "here\n";
  }

  FuncPtr get_solver(char filepath[])
  {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
    std::cout << "file path: " << filepath << "\n";
    return &solver;
  }

}  // namespace micm