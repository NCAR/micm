#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(BranchedConfig, DetectsInvalidConfig)
{
  micm::SolverConfig solver_config;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/branched/missing_reactants");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
  status = solver_config.ReadAndParse("./unit_configs/process/branched/missing_alkoxy_products");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
  status = solver_config.ReadAndParse("./unit_configs/process/branched/missing_nitrate_products");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
}

TEST(BranchedConfig, ParseConfig)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/branched/valid");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;
  EXPECT_EQ(process_vector.size(), 4);

  // first reaction
  {
    EXPECT_EQ(process_vector[0].reactants_.size(), 3);
    EXPECT_EQ(process_vector[0].reactants_[0].name_, "foo");
    EXPECT_EQ(process_vector[0].reactants_[1].name_, "quz");
    EXPECT_EQ(process_vector[0].reactants_[2].name_, "quz");
    EXPECT_EQ(process_vector[0].products_.size(), 2);
    EXPECT_EQ(process_vector[0].products_[0].first.name_, "bar");
    EXPECT_EQ(process_vector[0].products_[0].second, 1.0);
    EXPECT_EQ(process_vector[0].products_[1].first.name_, "baz");
    EXPECT_EQ(process_vector[0].products_[1].second, 3.2);
    micm::BranchedRateConstant* branched_rate_constant =
        dynamic_cast<micm::BranchedRateConstant*>(process_vector[0].rate_constant_.get());
    auto& params = branched_rate_constant->parameters_;
    EXPECT_EQ(params.X_, 12.3);
    EXPECT_EQ(params.Y_, 42.3);
    EXPECT_EQ(params.a0_, 1.0e-5);
    EXPECT_EQ(params.n_, 3);
    EXPECT_EQ(params.branch_, micm::BranchedRateConstantParameters::Branch::Alkoxy);
  } 
  {
    EXPECT_EQ(process_vector[1].reactants_.size(), 3);
    EXPECT_EQ(process_vector[1].reactants_[0].name_, "foo");
    EXPECT_EQ(process_vector[1].reactants_[1].name_, "quz");
    EXPECT_EQ(process_vector[1].reactants_[2].name_, "quz");
    EXPECT_EQ(process_vector[1].products_.size(), 1);
    EXPECT_EQ(process_vector[1].products_[0].first.name_, "quz");
    EXPECT_EQ(process_vector[1].products_[0].second, 1.0);
    micm::BranchedRateConstant* branched_rate_constant =
        dynamic_cast<micm::BranchedRateConstant*>(process_vector[1].rate_constant_.get());
    auto& params = branched_rate_constant->parameters_;
    EXPECT_EQ(params.X_, 12.3);
    EXPECT_EQ(params.Y_, 42.3);
    EXPECT_EQ(params.a0_, 1.0e-5);
    EXPECT_EQ(params.n_, 3);
    EXPECT_EQ(params.branch_, micm::BranchedRateConstantParameters::Branch::Nitrate);
  }

  // second reaction
  {
    EXPECT_EQ(process_vector[2].reactants_.size(), 2);
    EXPECT_EQ(process_vector[2].reactants_[0].name_, "bar");
    EXPECT_EQ(process_vector[2].reactants_[1].name_, "baz");
    EXPECT_EQ(process_vector[2].products_.size(), 1);
    EXPECT_EQ(process_vector[2].products_[0].first.name_, "baz");
    EXPECT_EQ(process_vector[2].products_[0].second, 1.0);
    micm::BranchedRateConstant* branched_rate_constant =
        dynamic_cast<micm::BranchedRateConstant*>(process_vector[2].rate_constant_.get());
    auto& params = branched_rate_constant->parameters_;
    EXPECT_EQ(params.X_, 0.32);
    EXPECT_EQ(params.Y_, 2.3e8);
    EXPECT_EQ(params.a0_, 0.423);
    EXPECT_EQ(params.n_, 6);
    EXPECT_EQ(params.branch_, micm::BranchedRateConstantParameters::Branch::Alkoxy);
  }
  {
    EXPECT_EQ(process_vector[3].reactants_.size(), 2);
    EXPECT_EQ(process_vector[3].reactants_[0].name_, "bar");
    EXPECT_EQ(process_vector[3].reactants_[1].name_, "baz");
    EXPECT_EQ(process_vector[3].products_.size(), 2);
    EXPECT_EQ(process_vector[3].products_[0].first.name_, "bar");
    EXPECT_EQ(process_vector[3].products_[0].second, 0.5);
    EXPECT_EQ(process_vector[3].products_[1].first.name_, "foo");
    EXPECT_EQ(process_vector[3].products_[1].second, 1.0);
    micm::BranchedRateConstant* branched_rate_constant =
        dynamic_cast<micm::BranchedRateConstant*>(process_vector[3].rate_constant_.get());
    auto& params = branched_rate_constant->parameters_;
    EXPECT_EQ(params.X_, 0.32);
    EXPECT_EQ(params.Y_, 2.3e8);
    EXPECT_EQ(params.a0_, 0.423);
    EXPECT_EQ(params.n_, 6);
    EXPECT_EQ(params.branch_, micm::BranchedRateConstantParameters::Branch::Nitrate);
  }
}
