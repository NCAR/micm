// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../dae_constraint_overshoot_policy.hpp"

#include <micm/Kokkos.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

using Builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>;

TEST(DAEConstraintOvershoot, AlgebraicVariableStaysNonNegative)
{
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  TestAlgebraicVariableStaysNonNegative(Builder(options));
}

TEST(DAEConstraintOvershoot, EquilibriumPlusConservation)
{
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  TestEquilibriumPlusConservation(Builder(options));
}

TEST(DAEConstraintOvershoot, AllRosenbrockOrdersConstrained)
{
  TestAllRosenbrockOrdersConstrained([](const micm::RosenbrockSolverParameters& options) { return Builder(options); });
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
