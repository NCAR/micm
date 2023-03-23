#ifdef USE_JSON
#include <micm/configure/system_builder.hpp>
#endif

#include <gtest/gtest.h>

#ifdef USE_JSON
TEST(SystemBuilder, DefaultConstructor){
  micm::SystemBuilder builder{};
}

TEST(SystemBuilder, DetectsInvalidConfigFile){
  micm::SystemBuilder builder{};
  EXPECT_ANY_THROW(builder.Build("not_a_config_file.json"));
}

TEST(SystemBuilder, DetectsInvalidConfigFileAndNoThrowDoesntThrow){
  micm::SystemBuilder<micm::JsonReaderPolicy, micm::NoThrowPolicy> builder{};
  EXPECT_NO_THROW(builder.Build("not_a_config_file.json"));

  auto system = builder.Build("not_a_config_file.json");
  ASSERT_TRUE(system.get() == nullptr);
}

TEST(SystemBuilder, JsonBuilder){
  micm::SystemBuilder builder{};
  auto system = builder.Build("unit_configs/chapman/config.json");

  EXPECT_TRUE(system != nullptr);
}
#endif