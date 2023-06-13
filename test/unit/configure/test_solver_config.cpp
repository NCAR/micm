#ifdef USE_JSON
#  include <micm/configure/solver_config.hpp>
#endif

#include <gtest/gtest.h>

#ifdef USE_JSON
TEST(SolverConfig, DetectsInvalidConfigFileAndThrow)
{
    micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig{};
    EXPECT_ANY_THROW(solverConfig.Configure("not_a_config_file.json"));
}

TEST(SolverConfig, DetectsInvalidConfigFileAndNoThrowDoesntThrow)
{
    micm::SolverConfig<micm::JsonReaderPolicy, micm::NoThrowPolicy> solverConfig{};
    EXPECT_NO_THROW(solverConfig.Configure("not_a_config_file.json"));

    std::variant<micm::SolverParameters, micm::ConfigErrorCode> configs = solverConfig.Configure("not_a_config_file.json");
    EXPECT_EQ(std::get<micm::ConfigErrorCode>(configs), micm::ConfigErrorCode::FileNotFound);
}

TEST(SolverConfig, ReadAndParse)
{
    micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig{};
    std::variant<micm::SolverParameters, micm::ConfigErrorCode> configs = solverConfig.Configure("./unit_configs/chapman/config.json");
    auto* solver_params_ptr = std::get_if<micm::SolverParameters>(&configs);

    EXPECT_TRUE(solver_params_ptr!= nullptr);
    
    micm::SolverParameters& solver_params = *solver_params_ptr;
    
    EXPECT_EQ(solver_params.system_.gas_phase_.species_.size(), 9);
    EXPECT_EQ(solver_params.processes_.size(), 7);

    int num_reactants_in_each_process[] = {1, 1, 1, 2, 2, 2, 3};
    short idx = 0;
    for (const auto& p : solver_params.processes_)
    {
        EXPECT_EQ(p.reactants_.size(), num_reactants_in_each_process[idx]);
        idx++;
    }
}
#endif