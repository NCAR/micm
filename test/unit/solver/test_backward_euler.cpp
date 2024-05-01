#include <micm/solver/backward_euler.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

TEST(BackwardEuler, DefaultConstructor)
{
  EXPECT_NO_THROW(micm::BackwardEuler be);
}
