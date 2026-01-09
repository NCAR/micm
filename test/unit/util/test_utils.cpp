// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/utils.hpp>

#include <gtest/gtest.h>

TEST(JoinStrings, EmptyVector)
{
  std::vector<std::string> names;
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "");
}

TEST(JoinStrings, SingleString)
{
  std::vector<std::string> names = { "species" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "species");
}

TEST(JoinStrings, TwoStrings)
{
  std::vector<std::string> names = { "accumulation", "aqueous" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "accumulation.aqueous");
}

TEST(JoinStrings, ThreeStrings)
{
  std::vector<std::string> names = { "accumulation", "aqueous", "CO2" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "accumulation.aqueous.CO2");
}

TEST(JoinStrings, EmptyStringAtStart)
{
  std::vector<std::string> names = { "", "accumulation", "aqueous" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "accumulation.aqueous");
}

TEST(JoinStrings, EmptyStringInMiddle)
{
  std::vector<std::string> names = { "accumulation", "", "aqueous" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "accumulation.aqueous");
}

TEST(JoinStrings, EmptyStringAtEnd)
{
  std::vector<std::string> names = { "accumulation", "aqueous", "" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "accumulation.aqueous");
}

TEST(JoinStrings, MultipleEmptyStrings)
{
  std::vector<std::string> names = { "", "accumulation", "", "aqueous", "" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "accumulation.aqueous");
}

TEST(JoinStrings, StringsWithDots)
{
  std::vector<std::string> names = { "accumulation.mode", "aqueous.phase" };
  std::string result = micm::JoinStrings(names);
  EXPECT_EQ(result, "accumulation.mode.aqueous.phase");
}

