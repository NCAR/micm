#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(SurfaceConfig, DetectsInvalidConfig)
{
  micm::SolverConfig solver_config;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/surface/missing_reactants");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
  status = solver_config.ReadAndParse("./unit_configs/process/surface/missing_products");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
  status = solver_config.ReadAndParse("./unit_configs/process/surface/missing_MUSICA_name");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
}

TEST(SurfaceConfig, ParseConfig)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/surface/valid");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);

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
    micm::SurfaceRateConstant* surface_rate_constant =
        dynamic_cast<micm::SurfaceRateConstant*>(process_vector[0].rate_constant_.get());
    EXPECT_EQ(surface_rate_constant->SizeCustomParameters(), 2);
    EXPECT_EQ(surface_rate_constant->CustomParameters()[0], "SURF.kfoo.effective radius [m]");
    EXPECT_EQ(surface_rate_constant->CustomParameters()[1], "SURF.kfoo.particle number concentration [# m-3]");
    EXPECT_EQ(surface_rate_constant->parameters_.reaction_probability_, 1.0);
    EXPECT_EQ(surface_rate_constant->diffusion_coefficient_, 2.3e-4);
    EXPECT_EQ(surface_rate_constant->mean_free_speed_factor_, 8.0 * GAS_CONSTANT / (M_PI * 0.123));
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
    micm::SurfaceRateConstant* surface_rate_constant =
        dynamic_cast<micm::SurfaceRateConstant*>(process_vector[1].rate_constant_.get());
    EXPECT_EQ(surface_rate_constant->SizeCustomParameters(), 2);
    EXPECT_EQ(surface_rate_constant->CustomParameters()[0], "SURF.kbar.effective radius [m]");
    EXPECT_EQ(surface_rate_constant->CustomParameters()[1], "SURF.kbar.particle number concentration [# m-3]");
    EXPECT_EQ(surface_rate_constant->parameters_.reaction_probability_, 0.5);
    EXPECT_EQ(surface_rate_constant->diffusion_coefficient_, 0.4e-5);
    EXPECT_EQ(surface_rate_constant->mean_free_speed_factor_, 8.0 * GAS_CONSTANT / (M_PI * 0.321));
  }
}

TEST(SurfaceConfig, DetectsNonstandardKeys)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/surface/contains_nonstandard_key");
  EXPECT_EQ(micm::ConfigParseStatus::ContainsNonStandardKey, status);
}

TEST(SurfaceConfig, DetectsNonstandardProductCoefficient)
{
  micm::SolverConfig solver_config;

  micm::ConfigParseStatus status = solver_config.ReadAndParse("./unit_configs/process/surface/nonstandard_product_coef");
  EXPECT_EQ(micm::ConfigParseStatus::ContainsNonStandardKey, status);
}