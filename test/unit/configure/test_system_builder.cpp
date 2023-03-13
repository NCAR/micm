#include <micm/configure/system_builder.hpp>

#include <gtest/gtest.h>

TEST(SystemBuilder, DefaultConstructor){
  micm::SystemBuilder builder{};
}

TEST(SystemBuilder, JsonBuilder){
  micm::SystemBuilder builder{};
  auto system = builder.Build("config/chapman.config");

  EXPECT_TRUE(system != nullptr);
}
