#include <micm/version.hpp>
#include <string>
#include <assert.h>

#include <gtest/gtest.h>

TEST(Version, FullVersion){
  auto version = getmicmVersion();
}

TEST(Version, VersionMajor){
  auto major = getmicmVersionMajor();
}

TEST(Version, VersionMinor){
  auto minor = getmicmVersionMinor();
}

TEST(Version, VersionPatch){
  auto patch = getmicmVersionPatch();
}

TEST(Version, VersionTweak){
 auto tweak = getmicmVersionTweak();
}