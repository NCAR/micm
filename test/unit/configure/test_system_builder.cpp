#include <micm/configure/system_builder.hpp>

#include <gtest/gtest.h>

TEST(SystemBuilder, DefaultConstructor){
  micm::SystemBuilder builder{};
}

TEST(SystemBuilder, DetectsInvalidConfigFile){
  micm::SystemBuilder builder{};
  EXPECT_ANY_THROW(builder.Build("not_a_config_file.json"));
}

TEST(SystemBuilder, JsonBuilder){
  micm::SystemBuilder builder{};
  auto system = builder.Build("unit_configs/chapman/config.json");

  EXPECT_TRUE(system != nullptr);
}
