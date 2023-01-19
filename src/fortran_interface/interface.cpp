#include <iostream>
#include "interface.hpp"

// namespace micm {

  void solver(double* arg1, double* arg2, double* result){ // NOLINT(misc-unused-parameters)
    std::cout << "here\n";
  }

  FuncPtr get_solver(char filepath[]){ // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
    std::cout << "file path: " << filepath << "\n";
    return &solver;
  }

// }