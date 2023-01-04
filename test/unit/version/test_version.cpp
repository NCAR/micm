#include <micm/version.hpp>

#include <gtest/gtest.h>

TEST(Version, FullVersion){
  auto version = micm::getmicmVersion();
}

TEST(Version, VersionMajor){
  auto major = micm::getmicmVersionMajor();
}

TEST(Version, VersionMinor){
  auto minor = micm::getmicmVersionMinor();
}

TEST(Version, VersionPatch){
  auto patch = micm::getmicmVersionPatch();
}

TEST(Version, VersionTweak){
 auto tweak = micm::getmicmVersionTweak();
}