/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#include <micm/micm.h>
#include <stdio.h>
#include <stdlib.h>

int main(const int argc, const char *argv[]) {

  if (argc != 2) {
    printf("\nUsage: c_host path/to/config.json\n\n");
    return 2;
  }

  micm_solver solver = create_solver(argv[1], 12);

  printf("\nNumber of cells: %d", solver.n_cells);
  printf("\nNumber of state variables: %d", solver.n_state_variables);
  for (int i=0; i<solver.n_state_variables; ++i) {
    printf("\n  species %d: '%s'", i, solver.state_variable_names[i]);
  }
  double *state = (double*) calloc(solver.n_state_variables, sizeof(double));
  printf("\nSolver return code: %d", solver.solve(state, 30.0));
  printf("\nAll done\n\n");

  return 0;
}
