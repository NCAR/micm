#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(System, DefaultConstructor){
  micm::System<double> system{};
}
