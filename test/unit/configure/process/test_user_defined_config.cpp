#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(UserDefinedConfig, DetectsInvalidConfig)
{
  micm::SolverConfig solver_config;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/user_defined/missing_reactants");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status) << " " << micm::configParseStatusToString(status) << "\n";
  status = solver_config.ReadAndParse("./unit_configs/process/user_defined/missing_products");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
  status = solver_config.ReadAndParse("./unit_configs/process/user_defined/missing_MUSICA_name");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
}

TEST(UserDefinedConfig, ParseConfig)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/user_defined/valid");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // first reaction
  {
    EXPECT_EQ(process_vector[0].reactants_.size(), 3);
    EXPECT_EQ(process_vector[0].reactants_[0].name_, "bar");
    EXPECT_EQ(process_vector[0].reactants_[1].name_, "bar");
    EXPECT_EQ(process_vector[0].reactants_[2].name_, "foo");
    EXPECT_EQ(process_vector[0].products_.size(), 2);
    EXPECT_EQ(process_vector[0].products_[0].first.name_, "baz");
    EXPECT_EQ(process_vector[0].products_[0].second, 1.4);
    EXPECT_EQ(process_vector[0].products_[1].first.name_, "foo");
    EXPECT_EQ(process_vector[0].products_[1].second, 1.0);
    micm::UserDefinedRateConstant* photo_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[0].rate_constant_.get());
    EXPECT_EQ(photo_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(photo_rate_constant->CustomParameters()[0], "USER.foo");
    EXPECT_EQ(photo_rate_constant->parameters_.scaling_factor_, 1.0);
  }

  // second reaction
  {
    EXPECT_EQ(process_vector[1].reactants_.size(), 2);
    EXPECT_EQ(process_vector[1].reactants_[0].name_, "foo");
    EXPECT_EQ(process_vector[1].reactants_[1].name_, "foo");
    EXPECT_EQ(process_vector[1].products_.size(), 1);
    EXPECT_EQ(process_vector[1].products_[0].first.name_, "bar");
    EXPECT_EQ(process_vector[1].products_[0].second, 1.0);
    micm::UserDefinedRateConstant* photo_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[1].rate_constant_.get());
    EXPECT_EQ(photo_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(photo_rate_constant->CustomParameters()[0], "USER.bar");
    EXPECT_EQ(photo_rate_constant->parameters_.scaling_factor_, 2.5);
  }
}

TEST(PhotolysisConfig, DetectsNonstandardKeys)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status =
      solver_config.ReadAndParse("./unit_configs/process/user_defined/contains_nonstandard_key");
  EXPECT_EQ(micm::ConfigParseStatus::ContainsNonStandardKey, status);
}

TEST(PhotolysisConfig, DetectsNonstandardProductCoefficient)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status =
      solver_config.ReadAndParse("./unit_configs/process/user_defined/nonstandard_product_coef");
  EXPECT_EQ(micm::ConfigParseStatus::ContainsNonStandardKey, status);
}

TEST(PhotolysisConfig, DetectsNonstandardReactantCoefficient)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status =
      solver_config.ReadAndParse("./unit_configs/process/user_defined/nonstandard_reactant_coef");
  EXPECT_EQ(micm::ConfigParseStatus::ContainsNonStandardKey, status);
}
