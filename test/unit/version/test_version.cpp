#include <gtest/gtest.h>

#include <micm/version.hpp>

TEST(Version, FullVersion)
{
  auto version = micm::getMicmVersion();
}

TEST(Version, VersionMajor)
{
  auto major = micm::getMicmVersionMajor();
}

TEST(Version, VersionMinor)
{
  auto minor = micm::getMicmVersionMinor();
}

TEST(Version, VersionPatch)
{
  auto patch = micm::getMicmVersionPatch();
}

TEST(Version, VersionTweak)
{
  auto tweak = micm::getMicmVersionTweak();
}