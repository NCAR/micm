// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "dae_algebraic_error_insensitivity_policy.hpp"

#include <micm/CPU.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

// The CPU driver keeps the standard (non-vectorized) matrix types this test has always used;
// the Kokkos driver uses the Kokkos vector-ordered ones. Parameterizing the shared header on
// the builder is what lets the two differ.
using Builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>;

namespace
{
  Builder MakeBuilder()
  {
    return Builder(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  }
}  // namespace

TEST(DAEAlgebraicError, ErrorSensitiveToBalanceAtol)
{
  TestErrorSensitiveToBalanceAtol([] { return MakeBuilder(); });
}

TEST(DAEAlgebraicError, AlgebraicVariableDoesNotOvershootDeeply)
{
  TestAlgebraicVariableDoesNotOvershootDeeply([] { return MakeBuilder(); });
}
