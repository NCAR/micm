/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_SOLVER_H
#define MICM_SOLVER_H

template <typename T> class Solver {
  private:
    const std::size_t size_;

  public:
    virtual void Solve(T[]) = 0;
};

#endif