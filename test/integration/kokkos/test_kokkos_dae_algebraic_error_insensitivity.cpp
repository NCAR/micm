// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../dae_algebraic_error_insensitivity_policy.hpp"

#include <micm/Kokkos.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

using Builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>;

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

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
