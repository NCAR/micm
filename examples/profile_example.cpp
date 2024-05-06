#include "util_example.cpp"

#include <micm/configure/solver_config.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using namespace micm;

template<template<class> class MatrixType, template<class> class SparseMatrixType>
int Run(const char* filepath, const char* initial_conditions, const std::string& matrix_ordering_type)
{
  using SolverType = RosenbrockSolver<MatrixType, SparseMatrixType>;

  fs::path config_path{ filepath };
  std::string initial_conditions_file{ initial_conditions };

  auto dataMap = readCSVToMap(initial_conditions_file);
  if (dataMap.empty())
  {
    std::cout << "Error: Failed in reading CSV file. Please verify file path and/or file name of \""
              << initial_conditions_file << "\"" << std::endl;
    return 1;
  }
  else if (dataMap.is_incomplete())
  {
    std::cout << "Error: Parsing csv file is incomplete . Please verify content formatting of \"" << initial_conditions_file
              << "\"" << std::endl;
    return 1;
  }

  SolverConfig solverConfig;
  solverConfig.ReadAndParse(config_path);

  SolverParameters solver_params = solverConfig.GetSolverParams();

  // add third-body species parameterizaton on air density
  for (auto& species : solver_params.system_.gas_phase_.species_)
    if (species.name_ == "M")
      species.parameterize_ = [](const Conditions& c) { return c.air_density_; };
  for (auto& process : solver_params.processes_)
  {
    for (auto& reactant : process.reactants_)
      if (reactant.name_ == "M")
        reactant.parameterize_ = [](const Conditions& c) { return c.air_density_; };
    for (auto& product : process.products_)
      if (product.first.name_ == "M")
        product.first.parameterize_ = [](const Conditions& c) { return c.air_density_; };
  }

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  auto params = RosenbrockSolverParameters::three_stage_rosenbrock_parameters();
  params.relative_tolerance_ = 0.1;
  RosenbrockSolver<> solver{ chemical_system, reactions, params };

  State state = solver.GetState();

  state.conditions_[0].temperature_ = dataMap.environments["temperature"];  // K
  state.conditions_[0].pressure_ = dataMap.environments["pressure"];        // Pa
  state.conditions_[0].air_density_ = dataMap.environments["air_density"];  // mol m-3

  state.SetConcentrations(dataMap.concentrations);

  if (!dataMap.custom_rate_params.empty())
  {
    for (auto& [rate_constant_key, rate_constant_val] : dataMap.custom_rate_params)
    {
      state.SetCustomRateParameter(rate_constant_key, rate_constant_val[0]);
    }
  }

  double time_step = dataMap.environments["time_step"];  // s
  double elapsed_solve_time = 0;
  typename SolverType::SolverResult result;

  MICM_PROFILE_BEGIN_SESSION("Runtime", "Profile-Runtime-" + matrix_ordering_type + ".json");
  while (elapsed_solve_time < time_step)
  {
    auto result = solver.Solve(time_step - elapsed_solve_time, state);
    elapsed_solve_time = result.final_time_;
    state.variables_ = result.result_;
    if (result.state_ != SolverState::Converged)
    {
      std::cout << "solver failed to converge" << std::endl;
      return 1;
    }
  }
  MICM_PROFILE_END_SESSION();

  return 0;
}

template<class T>
using SparseMatrixParam = micm::SparseMatrix<T>;
template<class T>
using Vector1MatrixParam = micm::VectorMatrix<T, 1>;
template<class T>
using Vector10MatrixParam = micm::VectorMatrix<T, 10>;
template<class T>
using Vector100MatrixParam = micm::VectorMatrix<T, 100>;
template<class T>
using Vector1000MatrixParam = micm::VectorMatrix<T, 1000>;
template<class T>
using Vector1SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Vector10SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<10>>;
template<class T>
using Vector100SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<100>>;
template<class T>
using Vector1000SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1000>>;

int main(const int argc, const char* argv[])
{
  if (argc < 3)
  {
    std::cout << "Usage: Profile MICM Solver <path> <initial condition>" << std::endl;
    std::cout << "  path                  Path to the folder holding the configuration data" << std::endl;
    std::cout << "  initial condition     csv file contanining initial condtions" << std::endl;
    return 1;
  }

  Run<Vector1000MatrixParam, Vector1000SparseMatrixParam>(argv[1], argv[2], "Vector-Sparse-1000");
};