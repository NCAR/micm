// Copyright (C) 2022 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
//
#include <micm/micm.h>
#include <micm/system/system.h>
#include <iostream>
#include <string.h>

#if __cplusplus
extern "C" {
#endif

unsigned int solver_function(double *state, double time_step__s) {
  std::cout << std::endl << "I'm solved!" << std::endl;
  return 0;
}

const micm_solver create_solver(const char* file_path, unsigned int n_cells) {
  std::cout << "file path: '" << file_path << "'" << std::endl;

  micm_solver solver{};
  solver.n_cells = n_cells;
  solver.n_state_variables = 3;
  solver.state_variable_names = (char**) malloc(sizeof(char*)*3);
  solver.state_variable_names[0] = (char*) malloc(sizeof(char)*2);
  strcpy(solver.state_variable_names[0],"O3");
  solver.state_variable_names[1] = (char*) malloc(sizeof(char)*2);
  strcpy(solver.state_variable_names[1], "NO");
  solver.state_variable_names[2] = (char*) malloc(sizeof(char)*3);
  strcpy(solver.state_variable_names[2], "NO2");
  solver.solve = solver_function;
  return solver;
}

#if __cplusplus
}
#endif
