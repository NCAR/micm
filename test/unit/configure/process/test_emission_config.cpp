#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(EmissionConfig, DetectsInvalidConfig)
{
  micm::SolverConfig solver_config;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/emission/missing_products");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
  status = solver_config.ReadAndParse("./unit_configs/process/emission/missing_MUSICA_name");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
}

TEST(EmissionConfig, ParseConfig)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/emission/valid");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // first reaction
  {
    EXPECT_EQ(process_vector[0].reactants_.size(), 0);
    EXPECT_EQ(process_vector[0].products_.size(), 1);
    EXPECT_EQ(process_vector[0].products_[0].first.name_, "foo");
    EXPECT_EQ(process_vector[0].products_[0].second, 1.0);
    micm::UserDefinedRateConstant* emission_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[0].rate_constant_.get());
    EXPECT_EQ(emission_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(emission_rate_constant->CustomParameters()[0], "EMIS.foo");
    EXPECT_EQ(emission_rate_constant->parameters_.scaling_factor_, 1.0);
  }

  // second reaction
  {
    EXPECT_EQ(process_vector[1].reactants_.size(), 0);
    EXPECT_EQ(process_vector[1].products_.size(), 1);
    EXPECT_EQ(process_vector[1].products_[0].first.name_, "bar");
    EXPECT_EQ(process_vector[1].products_[0].second, 1.0);
    micm::UserDefinedRateConstant* emission_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[1].rate_constant_.get());
    EXPECT_EQ(emission_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(emission_rate_constant->CustomParameters()[0], "EMIS.bar");
    EXPECT_EQ(emission_rate_constant->parameters_.scaling_factor_, 2.5);
  }
}

TEST(EmissionConfig, DetectsNonstandardKeys)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/emission/contains_nonstandard_key");
  EXPECT_EQ(micm::ConfigParseStatus::ContainsNonStandardKey, status);
}