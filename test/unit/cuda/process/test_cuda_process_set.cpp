#include <micm/cuda/process/cuda_process_set.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/process/process_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

using index_pair = std::pair<std::size_t, std::size_t>;

void compare_pair(const index_pair& a, const index_pair& b)
{
  EXPECT_EQ(a.first, b.first);
  EXPECT_EQ(a.second, b.second);
}

template<class CPUMatrixPolicy, class GPUMatrixPolicy>
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
  micm::State<CPUMatrixPolicy> cpu_state{ micm::StateParameters{
                                              .number_of_rate_constants_ = n_reactions,
                                              .variable_names_{ species_names },
                                              .custom_rate_parameter_labels_{},
                                          },
                                          n_cells };
  micm::State<GPUMatrixPolicy> gpu_state{ micm::StateParameters{
                                              .number_of_rate_constants_ = n_reactions,
                                              .variable_names_{ species_names },
                                              .custom_rate_parameter_labels_{},
                                          },
                                          n_cells };

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
    std::vector<micm::Yield> products{};
    for (std::size_t i_prod = 0; i_prod < n_product; ++i_prod)
    {
      products.push_back(micm::Yields(std::to_string(get_species_id()), 1.2));
    }
    processes.push_back(micm::Process::Create().SetReactants(reactants).SetProducts(products).SetPhase(gas_phase));
  }

  micm::ProcessSet cpu_set{ processes, cpu_state.variable_map_ };
  micm::CudaProcessSet gpu_set{ processes, gpu_state.variable_map_ };

  auto& cpu_state_vars = cpu_state.variables_.AsVector();
  std::ranges::generate(cpu_state_vars, get_double);
  gpu_state.variables_.AsVector().assign(cpu_state_vars.begin(), cpu_state_vars.end());
  gpu_state.variables_.CopyToDevice();

  CPUMatrixPolicy cpu_rate_constants{ n_cells, n_reactions };
  GPUMatrixPolicy gpu_rate_constants{ n_cells, n_reactions };
  auto& cpu_rate_vars = cpu_rate_constants.AsVector();
  std::ranges::generate(cpu_rate_vars, get_double);
  gpu_rate_constants.AsVector().assign(cpu_rate_vars.begin(), cpu_rate_vars.end());
  gpu_rate_constants.CopyToDevice();

  CPUMatrixPolicy cpu_forcing{ n_cells, n_species, 1000.0 };
  GPUMatrixPolicy gpu_forcing{ n_cells, n_species, 1000.0 };
  gpu_forcing.CopyToDevice();

  // kernel function call
  gpu_set.AddForcingTerms<GPUMatrixPolicy>(gpu_rate_constants, gpu_state.variables_, gpu_forcing);
  gpu_forcing.CopyToHost();

  // CPU function call
  cpu_set.AddForcingTerms<CPUMatrixPolicy>(cpu_rate_constants, cpu_state.variables_, cpu_forcing);

  // checking accuracy with comparison between CPU and GPU result
  std::vector<double> cpu_forcing_vector = cpu_forcing.AsVector();
  std::vector<double> gpu_forcing_vector = gpu_forcing.AsVector();

  for (int i = 0; i < cpu_forcing_vector.size(); i++)
  {
    double a = cpu_forcing_vector[i];
    double b = gpu_forcing_vector[i];
    EXPECT_LT(std::abs((a - b) / a), 1.e-11);
  }
}

template<class CPUMatrixPolicy, class CPUSparseMatrixPolicy, class GPUDenseMatrixPolicy, class GPUSparseMatrixPolicy>
void testRandomSystemSubtractJacobianTerms(std::size_t n_cells, std::size_t n_reactions, std::size_t n_species)
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

  micm::State<CPUMatrixPolicy> cpu_state{ micm::StateParameters{
                                              .number_of_rate_constants_ = n_reactions,
                                              .variable_names_{ species_names },
                                              .custom_rate_parameter_labels_{},
                                          },
                                          n_cells };
  micm::State<GPUDenseMatrixPolicy> gpu_state{ micm::StateParameters{
                                                   .number_of_rate_constants_ = n_reactions,
                                                   .variable_names_{ species_names },
                                                   .custom_rate_parameter_labels_{},
                                               },
                                               n_cells };

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
    std::vector<micm::Yield> products{};
    for (std::size_t i_prod = 0; i_prod < n_product; ++i_prod)
    {
      products.push_back(micm::Yields(std::to_string(get_species_id()), 1.2));
    }
    processes.push_back(micm::Process::Create().SetReactants(reactants).SetProducts(products).SetPhase(gas_phase));
  }

  micm::ProcessSet cpu_set{ processes, cpu_state.variable_map_ };
  micm::CudaProcessSet gpu_set{ processes, gpu_state.variable_map_ };

  auto& cpu_state_vars = cpu_state.variables_.AsVector();
  std::ranges::generate(cpu_state_vars, get_double);
  gpu_state.variables_.AsVector().assign(cpu_state_vars.begin(), cpu_state_vars.end());
  gpu_state.variables_.CopyToDevice();

  CPUMatrixPolicy cpu_rate_constants{ n_cells, n_reactions };
  GPUDenseMatrixPolicy gpu_rate_constants{ n_cells, n_reactions };
  auto& cpu_rate_vars = cpu_rate_constants.AsVector();
  std::ranges::generate(cpu_rate_vars, get_double);
  gpu_rate_constants.AsVector().assign(cpu_rate_vars.begin(), cpu_rate_vars.end());
  gpu_rate_constants.CopyToDevice();

  auto non_zero_elements = cpu_set.NonZeroJacobianElements();

  auto cpu_builder = CPUSparseMatrixPolicy::Create(n_species).SetNumberOfBlocks(n_cells).InitialValue(100.0);
  for (auto& elem : non_zero_elements)
    cpu_builder = cpu_builder.WithElement(elem.first, elem.second);
  CPUSparseMatrixPolicy cpu_jacobian{ cpu_builder };
  auto gpu_builder = GPUSparseMatrixPolicy::Create(n_species).SetNumberOfBlocks(n_cells).InitialValue(100.0);
  for (auto& elem : non_zero_elements)
    gpu_builder = gpu_builder.WithElement(elem.first, elem.second);
  GPUSparseMatrixPolicy gpu_jacobian{ gpu_builder };
  gpu_jacobian.CopyToDevice();

  cpu_set.SetJacobianFlatIds(cpu_jacobian);
  gpu_set.SetJacobianFlatIds(gpu_jacobian);

  cpu_set.SubtractJacobianTerms(cpu_rate_constants, cpu_state.variables_, cpu_jacobian);
  gpu_set.SubtractJacobianTerms(gpu_rate_constants, gpu_state.variables_, gpu_jacobian);
  gpu_jacobian.CopyToHost();

  // checking accuracy of jacobian between CPU and GPU
  std::vector<double> cpu_jacobian_vector = cpu_jacobian.AsVector();
  std::vector<double> gpu_jacobian_vector = gpu_jacobian.AsVector();

  for (int i = 0; i < cpu_jacobian_vector.size(); i++)
  {
    double a = cpu_jacobian_vector[i];
    double b = gpu_jacobian_vector[i];
    EXPECT_LT(std::abs((a - b) / a), 2.e-10);
  }
}

using Group10000VectorMatrix = micm::VectorMatrix<double, 10000>;
using Group10000SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<10000>>;
using Group10000CudaDenseMatrix = micm::CudaDenseMatrix<double, 10000>;
using Group10000CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<10000>>;

TEST(RandomCudaProcessSet, Forcing)
{
  testRandomSystemAddForcingTerms<Group10000VectorMatrix, Group10000CudaDenseMatrix>(10000, 500, 400);
}
TEST(RandomCudaProcessSet, Jacobian)
{
  testRandomSystemSubtractJacobianTerms<
      Group10000VectorMatrix,
      Group10000SparseVectorMatrix,
      Group10000CudaDenseMatrix,
      Group10000CudaSparseMatrix>(10000, 500, 400);
}
