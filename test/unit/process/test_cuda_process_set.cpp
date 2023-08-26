#include <gtest/gtest.h>

#include <chrono>
#include <functional>
#include <iostream>
#include <micm/process/cuda_process_set.hpp>
#include <micm/process/process_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>
#include <vector>

using yields = std::pair<micm::Species, double>;
using index_pair = std::pair<std::size_t, std::size_t>;

void compare_pair(const index_pair& a, const index_pair& b)
{
  EXPECT_EQ(a.first, b.first);
  EXPECT_EQ(a.second, b.second);
}

template<template<class> class MatrixPolicy>
void testRandomSystemAddForcingTerms(std::size_t n_cells, std::size_t n_reactions, std::size_t n_species)
{
  auto get_n_react = std::bind(std::uniform_int_distribution<>(0, 3), std::default_random_engine());
  auto get_n_product = std::bind(std::uniform_int_distribution<>(0, 10), std::default_random_engine());
  auto get_species_id = std::bind(std::uniform_int_distribution<>(0, n_species - 1), std::default_random_engine());
  std::lognormal_distribution<double> distribution(-2.0, 4.0);
  auto get_double = std::bind(distribution, std::default_random_engine());

  std::vector<micm::Species> species{};
  std::vector<std::string> species_names{};
  for (std::size_t i = 0; i < n_species; ++i)
  {
    species.push_back(micm::Species{ std::to_string(i) });
    species_names.push_back(std::to_string(i));
  }
  micm::Phase gas_phase{ species };
  micm::State<MatrixPolicy> state{ micm::StateParameters{ .state_variable_names_{ species_names },
                                                          .custom_rate_parameter_labels_{},
                                                          .number_of_grid_cells_ = n_cells,
                                                          .number_of_rate_constants_ = n_reactions } };

  std::vector<micm::Process> processes{};
  for (std::size_t i = 0; i < n_reactions; ++i)
  {
    auto n_react = get_n_react();
    std::vector<micm::Species> reactants{};
    for (std::size_t i_react = 0; i_react < n_react; ++i_react)
    {
      reactants.push_back({ std::to_string(get_species_id()) });
    }
    auto n_product = get_n_product();
    std::vector<yields> products{};
    for (std::size_t i_prod = 0; i_prod < n_product; ++i_prod)
    {
      products.push_back(yields(std::to_string(get_species_id()), 1.2));
    }
    processes.push_back(micm::Process::create().reactants(reactants).products(products).phase(gas_phase));
  }

  micm::ProcessSet cpu_set{ processes, state };
  micm::CudaProcessSet gpu_set{ processes, state };

  for (auto& elem : state.variables_.AsVector())
    elem = get_double();

  MatrixPolicy<double> rate_constants{ n_cells, n_reactions };
  for (auto& elem : rate_constants.AsVector())
    elem = get_double();

  MatrixPolicy<double> cpu_forcing{ n_cells, n_species, 1000.0 };
  MatrixPolicy<double> gpu_forcing{};
  gpu_forcing = cpu_forcing;

  // kernel function call
  gpu_set.AddForcingTerms<MatrixPolicy>(rate_constants, state.variables_, gpu_forcing);

  // CPU function call
  cpu_set.AddForcingTerms<MatrixPolicy>(rate_constants, state.variables_, cpu_forcing);

  // checking accuracy with comparison between CPU and GPU result
  std::vector<double> cpu_forcing_vector = cpu_forcing.AsVector();
  std::vector<double> gpu_forcing_vector = gpu_forcing.AsVector();

  for (int i = 0; i < cpu_forcing_vector.size(); i++)
  {
    double a = cpu_forcing_vector[i];
    double b = gpu_forcing_vector[i];
    EXPECT_EQ(a, b);
  }
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void testRandomSystemAddJacobianTerms(std::size_t n_cells, std::size_t n_reactions, std::size_t n_species)
{
  auto get_n_react = std::bind(std::uniform_int_distribution<>(0, 3), std::default_random_engine());
  auto get_n_product = std::bind(std::uniform_int_distribution<>(0, 10), std::default_random_engine());
  auto get_species_id = std::bind(std::uniform_int_distribution<>(0, n_species - 1), std::default_random_engine());
  std::lognormal_distribution<double> distribution(-2.0, 4.0);
  auto get_double = std::bind(distribution, std::default_random_engine());

  std::vector<micm::Species> species{};
  std::vector<std::string> species_names{};
  for (std::size_t i = 0; i < n_species; ++i)
  {
    species.push_back(micm::Species{ std::to_string(i) });
    species_names.push_back(std::to_string(i));
  }
  micm::Phase gas_phase{ species };
  micm::State<MatrixPolicy> state{ micm::StateParameters{ .state_variable_names_{ species_names },
                                                          .number_of_grid_cells_ = n_cells,    
                                                          .number_of_rate_constants_ = n_reactions } };

  std::vector<micm::Process> processes{};
  for (std::size_t i = 0; i < n_reactions; ++i)
  {
    auto n_react = get_n_react();
    std::vector<micm::Species> reactants{};
    for (std::size_t i_react = 0; i_react < n_react; ++i_react)
    {
      reactants.push_back({ std::to_string(get_species_id()) });
    }
    auto n_product = get_n_product();
    std::vector<yields> products{};
    for (std::size_t i_prod = 0; i_prod < n_product; ++i_prod)
    {
      products.push_back(yields(std::to_string(get_species_id()), 1.2));
    }
    processes.push_back(micm::Process::create().reactants(reactants).products(products).phase(gas_phase));
  }

  micm::ProcessSet cpu_set{ processes, state };
  micm::CudaProcessSet gpu_set{ processes, state };

  for (auto& elem : state.variables_.AsVector())
    elem = get_double();

  MatrixPolicy<double> rate_constants{ n_cells, n_reactions };
  for (auto& elem : rate_constants.AsVector())
    elem = get_double();

  auto non_zero_elements = cpu_set.NonZeroJacobianElements();
  auto builder = SparseMatrixPolicy<double>::create(n_species).number_of_blocks(n_cells).initial_value(100.0);
  for (auto& elem : non_zero_elements)
    builder = builder.with_element(elem.first, elem.second);
  SparseMatrixPolicy<double> cpu_jacobian{ builder };
  SparseMatrixPolicy<double> gpu_jacobian{ builder };

  cpu_set.SetJacobianFlatIds(cpu_jacobian);
  gpu_set.SetJacobianFlatIds(gpu_jacobian);

  cpu_set.AddJacobianTerms<MatrixPolicy, SparseMatrixPolicy>(rate_constants, state.variables_, cpu_jacobian);
  gpu_set.AddJacobianTerms<MatrixPolicy, SparseMatrixPolicy>(rate_constants, state.variables_, gpu_jacobian);

  // checking accuracy of jacobian between CPU and GPU
  std::vector<double> cpu_jacobian_vector = cpu_jacobian.AsVector();
  std::vector<double> gpu_jacobian_vector = gpu_jacobian.AsVector();

  for (int i = 0; i < cpu_jacobian_vector.size(); i++)
  {
    double a = cpu_jacobian_vector[i];
    double b = gpu_jacobian_vector[i];
    ASSERT_EQ(a, b);
  }
}

template<class T>
using Group10000VectorMatrix = micm::VectorMatrix<T, 5>;

template<class T>
using Group10000SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<5>>;

TEST(RandomCudaProcessSet, Forcing)
{
  testRandomSystemAddForcingTerms<Group10000VectorMatrix>(5, 8, 6);
}
TEST(RandomCudaProcessSet, Jacobian)
{
  testRandomSystemAddJacobianTerms<Group10000VectorMatrix, Group10000SparseVectorMatrix>(5, 8,6);
}
