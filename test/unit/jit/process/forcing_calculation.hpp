// A function hard-coded to calculate forcing for a toy chemistry system.
// This can be used to compare with the corresponding JITed function

#pragma once

#define NUM_GRID_CELLS 1000

void calculate_forcing(double* rate_constants, double* state, double* forcing);