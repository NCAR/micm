#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(PhotolysisConfig, DetectsInvalidConfig)
{
  micm::SolverConfig solver_config;

  // Read and parse the configure files
  try {
    solver_config.ReadAndParse("./unit_configs/process/photolysis/missing_reactants");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::RequiredKeyNotFound));
  }
  try {
    solver_config.ReadAndParse("./unit_configs/process/photolysis/missing_products");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::RequiredKeyNotFound));
  }
  try {
    solver_config.ReadAndParse("./unit_configs/process/photolysis/missing_MUSICA_name");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::RequiredKeyNotFound));
  }
}

TEST(PhotolysisConfig, ParseConfig)
{
  micm::SolverConfig solver_config;

  EXPECT_NO_THROW(solver_config.ReadAndParse("./unit_configs/process/photolysis/valid"));

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // first reaction
  {
    EXPECT_EQ(process_vector[0].reactants_.size(), 1);
    EXPECT_EQ(process_vector[0].reactants_[0].name_, "foo");
    EXPECT_EQ(process_vector[0].products_.size(), 2);
    EXPECT_EQ(process_vector[0].products_[0].first.name_, "bar");
    EXPECT_EQ(process_vector[0].products_[0].second, 1.0);
    EXPECT_EQ(process_vector[0].products_[1].first.name_, "baz");
    EXPECT_EQ(process_vector[0].products_[1].second, 3.2);
    micm::UserDefinedRateConstant* photo_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[0].rate_constant_.get());
    EXPECT_EQ(photo_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(photo_rate_constant->CustomParameters()[0], "PHOTO.jfoo");
    EXPECT_EQ(photo_rate_constant->parameters_.scaling_factor_, 1.0);
  }

  // second reaction
  {
    EXPECT_EQ(process_vector[1].reactants_.size(), 1);
    EXPECT_EQ(process_vector[1].reactants_[0].name_, "bar");
    EXPECT_EQ(process_vector[1].products_.size(), 2);
    EXPECT_EQ(process_vector[1].products_[0].first.name_, "bar");
    EXPECT_EQ(process_vector[1].products_[0].second, 0.5);
    EXPECT_EQ(process_vector[1].products_[1].first.name_, "foo");
    EXPECT_EQ(process_vector[1].products_[1].second, 1.0);
    micm::UserDefinedRateConstant* photo_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[1].rate_constant_.get());
    EXPECT_EQ(photo_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(photo_rate_constant->CustomParameters()[0], "PHOTO.jbar");
    EXPECT_EQ(photo_rate_constant->parameters_.scaling_factor_, 2.5);
  }
}

TEST(PhotolysisConfig, DetectsNonstandardKeys)
{
  micm::SolverConfig solver_config;

  try {
    solver_config.ReadAndParse("./unit_configs/process/photolysis/contains_nonstandard_key");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::ContainsNonStandardKey));
  }
}

TEST(PhotolysisConfig, DetectsNonstandardProductCoefficient)
{
  micm::SolverConfig solver_config;

  try {
    solver_config.ReadAndParse("./unit_configs/process/photolysis/nonstandard_product_coef");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::ContainsNonStandardKey));
  }
}

TEST(PhotolysisConfig, DetectsNonstandardReactantCoefficient)
{
  micm::SolverConfig solver_config;

  try {
    solver_config.ReadAndParse("./unit_configs/process/photolysis/nonstandard_reactant_coef");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::ContainsNonStandardKey));
  }
}