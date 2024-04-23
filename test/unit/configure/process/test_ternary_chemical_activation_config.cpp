#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(TernaryChemicalActivationConfig, DetectsInvalidConfig)
{
  // Read and parse the configure files
  micm::SolverConfig solver_config;
  try {
    solver_config.ReadAndParse("./unit_configs/process/ternary_chemical_activation/missing_reactants");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::RequiredKeyNotFound));
  }
  try {
    solver_config.ReadAndParse("./unit_configs/process/ternary_chemical_activation/missing_products");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::RequiredKeyNotFound));
  }
}

TEST(TernaryChemicalActivationConfig, ParseConfig)
{
  micm::SolverConfig solver_config;
  EXPECT_NO_THROW(solver_config.ReadAndParse("./unit_configs/process/ternary_chemical_activation/valid"));

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // Convert Arrhenius parameters from expecting molecules cm-3 to moles m-3
  const double conv = 1.0e-6 * 6.02214076e23;

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
    micm::TernaryChemicalActivationRateConstant* ternary_rate_constant =
        dynamic_cast<micm::TernaryChemicalActivationRateConstant*>(process_vector[0].rate_constant_.get());
    auto& params = ternary_rate_constant->parameters_;
    EXPECT_EQ(params.k0_A_, 1.0 * conv * conv);
    EXPECT_EQ(params.k0_B_, 0.0);
    EXPECT_EQ(params.k0_C_, 0.0);
    EXPECT_EQ(params.kinf_A_, 1.0 * conv);
    EXPECT_EQ(params.kinf_B_, 0.0);
    EXPECT_EQ(params.kinf_C_, 0.0);
    EXPECT_EQ(params.Fc_, 0.6);
    EXPECT_EQ(params.N_, 1.0);
  }

  // second reaction
  {
    EXPECT_EQ(process_vector[1].reactants_.size(), 2);
    EXPECT_EQ(process_vector[1].reactants_[0].name_, "bar");
    EXPECT_EQ(process_vector[1].reactants_[1].name_, "baz");
    EXPECT_EQ(process_vector[1].products_.size(), 2);
    EXPECT_EQ(process_vector[1].products_[0].first.name_, "bar");
    EXPECT_EQ(process_vector[1].products_[0].second, 0.5);
    EXPECT_EQ(process_vector[1].products_[1].first.name_, "foo");
    EXPECT_EQ(process_vector[1].products_[1].second, 1.0);
    micm::TernaryChemicalActivationRateConstant* ternary_rate_constant =
        dynamic_cast<micm::TernaryChemicalActivationRateConstant*>(process_vector[1].rate_constant_.get());
    auto& params = ternary_rate_constant->parameters_;
    EXPECT_EQ(params.k0_A_, 32.1 * conv);
    EXPECT_EQ(params.k0_B_, -2.3);
    EXPECT_EQ(params.k0_C_, 102.3);
    EXPECT_EQ(params.kinf_A_, 63.4);
    EXPECT_EQ(params.kinf_B_, -1.3);
    EXPECT_EQ(params.kinf_C_, 908.5);
    EXPECT_EQ(params.Fc_, 1.3);
    EXPECT_EQ(params.N_, 32.1);
  }
}

TEST(TernaryChemicalActivationConfig, DetectsNonstandardKeys)
{
  micm::SolverConfig solver_config;
  try {
    solver_config.ReadAndParse("./unit_configs/process/ternary_chemical_activation/contains_nonstandard_key");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::ContainsNonStandardKey));
  }
}

TEST(TernaryChemicalActivationConfig, DetectsNonstandardProductCoefficient)
{
  micm::SolverConfig solver_config;
  try {
    solver_config.ReadAndParse("./unit_configs/process/ternary_chemical_activation/nonstandard_product_coef");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::ContainsNonStandardKey));
  }
}

TEST(TernaryChemicalActivationConfig, DetectsNonstandardReactantCoefficient)
{
  micm::SolverConfig solver_config;
  try {
    solver_config.ReadAndParse("./unit_configs/process/ternary_chemical_activation/nonstandard_reactant_coef");
  } catch (const std::system_error& e) {
    EXPECT_EQ(e.code().value(), static_cast<int>(MicmConfigErrc::ContainsNonStandardKey));
  }
}