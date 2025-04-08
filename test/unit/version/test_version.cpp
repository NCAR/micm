#include <micm/version.hpp>

#include <gtest/gtest.h>

TEST(Version, FullVersion)
{
  auto version = micm::GetMicmVersion();
}

TEST(Version, VersionMajor)
{
  auto major = micm::GetMicmVersionMajor();
}

TEST(Version, VersionMinor)
{
  auto minor = micm::GetMicmVersionMinor();
}

TEST(Version, VersionPatch)
{
  auto patch = micm::GetMicmVersionPatch();
}

TEST(Version, VersionTweak)
{
  auto tweak = micm::GetMicmVersionTweak();
}