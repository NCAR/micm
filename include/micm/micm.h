/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_MICM_H
#define MICM_MICM_H

#ifdef __cplusplus
extern "C" {
#endif

/// Solver and related functions that can be used by a host-model to
/// advance the chemical system in time
typedef struct micm_solver {
  /// Number of grid cells to be solved
  int n_cells;
  /// Number of state variables per grid cell
  int n_state_variables;
  /// Unique names for each state variable per grid cell
  char **state_variable_names;
  /// Pointer to the solver function
  ///
  /// state holds the solver state variables in the order specified in
  /// micm_solver->state_variable_names[] for each independent grid
  /// cell requested in the creation of the solver
  ///
  /// time_step__s is the time in seconds to advance the state
  ///
  /// The return value is a code where 0 indicates success and
  /// any other value indicates a solver failure
  unsigned int (*solve)(double *state, double time_step__s);
} micm_solver;

/// Returns a function that can be used to advance the chemical system
///
/// file_path is the path to the configuration data for the chemical
///           system
/// n_cells is the number of grid cells that will be solved simultaneously
///         by the generated solver
///
/// Returns a micm_solver object created for the specified system
const micm_solver create_solver(const char *file_path, unsigned int n_cells);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
