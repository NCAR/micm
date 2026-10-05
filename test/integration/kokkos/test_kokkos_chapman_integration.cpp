// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../chapman_policy.hpp"

#include <micm/Kokkos.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

TEST(ChapmanIntegration, CanBuildChapmanSystem)
{
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  TestChapmanSystem(micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(options));
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
