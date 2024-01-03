// Copyright (C) 2023-2024 National Center for Atmospheric Research,
// SPDX-License-Identifier: Apache-2.0
//
// A function hard-coded to calculate forcing for a toy chemistry system.
// This can be used to compare with the corresponding JITed function

#include "forcing_calculation.hpp"

#define SPEC_A 0
#define SPEC_B 1
#define SPEC_C 2
#define SPEC_D 3
#define SPEC_E 4
#define SPEC_F 5

void calculate_forcing(double* rate_constants, double* state, double* forcing)
{
  double rate[NUM_GRID_CELLS];

  // Reaction 1: A + B -> 3.2 D
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] = rate_constants[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] *= state[SPEC_A * NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] *= state[SPEC_B * NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_A * NUM_GRID_CELLS + i] -= rate[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_B * NUM_GRID_CELLS + i] -= rate[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_D * NUM_GRID_CELLS + i] += 3.2 * rate[i];

  // Reaction 2: E + C -> A + 2 F
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] = rate_constants[NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] *= state[SPEC_E * NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] *= state[SPEC_C * NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_E * NUM_GRID_CELLS + i] -= rate[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_C * NUM_GRID_CELLS + i] -= rate[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_A * NUM_GRID_CELLS + i] += rate[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_F * NUM_GRID_CELLS + i] += 2.0 * rate[i];

  // Reaction 3: C + B -> A
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] = rate_constants[2 * NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] *= state[SPEC_C * NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    rate[i] *= state[SPEC_B * NUM_GRID_CELLS + i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_C * NUM_GRID_CELLS + i] -= rate[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_B * NUM_GRID_CELLS + i] -= rate[i];
  for (int i = 0; i < NUM_GRID_CELLS; ++i)
    forcing[SPEC_A * NUM_GRID_CELLS + i] += rate[i];
}