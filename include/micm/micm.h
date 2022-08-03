/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_MICM_H
#define MICM_MICM_H

#ifdef __cplusplus
extern "C" {
#endif

/// Returns a function that can be used to advance the chemical system
void create_solver(char *file_path);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
