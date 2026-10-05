// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "chapman_policy.hpp"

#include <micm/CPU.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

TEST(ChapmanIntegration, CanBuildChapmanSystem)
{
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  TestChapmanSystem(micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options));
}
