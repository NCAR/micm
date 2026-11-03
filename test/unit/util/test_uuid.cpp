// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/uuid.hpp>

#include <gtest/gtest.h>

#include <regex>
#include <string>
#include <unordered_set>

TEST(Uuid, MatchesVersion4Format)
{
  const std::regex uuid_v4("^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$");
  for (int i = 0; i < 100; ++i)
  {
    auto uuid = micm::GenerateUuid();
    EXPECT_EQ(uuid.size(), 36);
    EXPECT_TRUE(std::regex_match(uuid, uuid_v4)) << uuid;
  }
}

TEST(Uuid, GeneratesUniqueValues)
{
  std::unordered_set<std::string> uuids;
  const int n = 10000;
  for (int i = 0; i < n; ++i)
  {
    uuids.insert(micm::GenerateUuid());
  }
  EXPECT_EQ(uuids.size(), n);
}
