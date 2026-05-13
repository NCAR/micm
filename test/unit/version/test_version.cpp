#include <micm/version.hpp>

#include <iostream>

#include <gtest/gtest.h>

TEST(Version, FullVersion)
{
  auto version = micm::GetMicmVersion();
  std::cout << "Micm version: " << version << std::endl;
}

TEST(Version, VersionMajor)
{
  auto major = micm::GetMicmVersionMajor();
  std::cout << "Micm version major: " << major << std::endl;
}

TEST(Version, VersionMinor)
{
  auto minor = micm::GetMicmVersionMinor();
  std::cout << "Micm version minor: " << minor << std::endl;
}

TEST(Version, VersionPatch)
{
  auto patch = micm::GetMicmVersionPatch();
  std::cout << "Micm version patch: " << patch << std::endl;
}

TEST(Version, VersionTweak)
{
  auto tweak = micm::GetMicmVersionTweak();
  std::cout << "Micm version tweak: " << tweak << std::endl;
}