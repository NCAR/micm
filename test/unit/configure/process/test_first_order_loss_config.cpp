#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(FirstOrderLossConfig, DetectsInvalidConfig)
{
  micm::SolverConfig solver_config;

  // Read and parse the configure files
  try
  {
    solver_config.ReadAndParse("./unit_configs/process/first_order_loss/missing_reactants");
  }
  catch (const std::system_error& e)
  {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::RequiredKeyNotFound));
  }
  try
  {
    solver_config.ReadAndParse("./unit_configs/process/first_order_loss/missing_MUSICA_name");
  }
  catch (const std::system_error& e)
  {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::RequiredKeyNotFound));
  }
}

TEST(FirstOrderLossConfig, ParseConfig)
{
  micm::SolverConfig solver_config;

  EXPECT_NO_THROW(solver_config.ReadAndParse("./unit_configs/process/first_order_loss/valid"));

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // first reaction
  {
    EXPECT_EQ(process_vector[0].reactants_.size(), 1);
    EXPECT_EQ(process_vector[0].reactants_[0].name_, "foo");
    EXPECT_EQ(process_vector[0].products_.size(), 0);
    micm::UserDefinedRateConstant* first_order_loss_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[0].rate_constant_.get());
    EXPECT_EQ(first_order_loss_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(first_order_loss_rate_constant->CustomParameters()[0], "LOSS.foo");
    EXPECT_EQ(first_order_loss_rate_constant->parameters_.scaling_factor_, 1.0);
  }

  // second reaction
  {
    EXPECT_EQ(process_vector[1].reactants_.size(), 1);
    EXPECT_EQ(process_vector[1].reactants_[0].name_, "bar");
    EXPECT_EQ(process_vector[1].products_.size(), 0);
    micm::UserDefinedRateConstant* first_order_loss_rate_constant =
        dynamic_cast<micm::UserDefinedRateConstant*>(process_vector[1].rate_constant_.get());
    EXPECT_EQ(first_order_loss_rate_constant->SizeCustomParameters(), 1);
    EXPECT_EQ(first_order_loss_rate_constant->CustomParameters()[0], "LOSS.bar");
    EXPECT_EQ(first_order_loss_rate_constant->parameters_.scaling_factor_, 2.5);
  }
}

TEST(FirstOrderLossConfig, DetectsNonstandardKeys)
{
  micm::SolverConfig solver_config;

  try
  {
    solver_config.ReadAndParse("./unit_configs/process/first_order_loss/contains_nonstandard_key");
  }
  catch (const std::system_error& e)
  {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::ContainsNonStandardKey));
  }
}