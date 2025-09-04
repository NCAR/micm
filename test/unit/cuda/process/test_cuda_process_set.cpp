#include <micm/cuda/process/cuda_process_set.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/process_set.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
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
  micm::ArrheniusRateConstant arrhenius_rate_constant({ .A_ = 12.2, .C_ = 300.0 });

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
      products.push_back(micm::Yield(std::to_string(get_species_id()), 1.2));
    }
    processes.push_back(micm::ChemicalReactionBuilder()
                            .SetReactants(reactants)
                            .SetProducts(products)
                            .SetPhase(gas_phase)
                            .SetRateConstant(arrhenius_rate_constant)
                            .Build());
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

  for (int i = 0; i < n_cells; i++)
    for (int j = 0; j < n_species; j++)
    {
      double a = cpu_forcing[i][j];
      double b = gpu_forcing[i][j];
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
  micm::ArrheniusRateConstant arrhenius_rate_constant({ .A_ = 12.2, .C_ = 300.0 });

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
      products.push_back(micm::Yield(std::to_string(get_species_id()), 1.2));
    }
    processes.push_back(micm::ChemicalReactionBuilder()
                            .SetReactants(reactants)
                            .SetProducts(products)
                            .SetRateConstant(arrhenius_rate_constant)
                            .SetPhase(gas_phase)
                            .Build());
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


  for (int i = 0; i < n_cells; i++)
    for (int j = 0; j < n_species; j++)
      for (int k = 0; k < n_species; k++)
      {
        if (cpu_jacobian.IsZero(j, k))
          continue;
        ASSERT_TRUE(!gpu_jacobian.IsZero(j, k));
        double a = cpu_jacobian[i][j][k];
        double b = gpu_jacobian[i][j][k];
        EXPECT_LT(std::abs((a - b) / a), 2.e-10);
      }
}

using FloatingPointType = double;

using Group1CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 1>;
using Group3CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 3>;
using Group27CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 27>;
using Group32CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 32>;
using Group43CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 43>;
using Group77CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 77>;
using Group113CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 113>;
using Group193CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 193>;
using Group281CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 281>;
using Group472CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 472>;
using Group512CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 512>;
using Group739CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 739>;
using Group1130CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 1130>;

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group3CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group27CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<27>>;
using Group32CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<32>>;
using Group43CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<43>>;
using Group77CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<77>>;
using Group113CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<113>>;
using Group193CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<193>>;
using Group281CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<281>>;
using Group472CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<472>>;
using Group512CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<512>>;
using Group739CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<739>>;
using Group1130CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1130>>;

using Group1VectorMatrix = micm::VectorMatrix<FloatingPointType, 1>;
using Group3VectorMatrix = micm::VectorMatrix<FloatingPointType, 3>;
using Group27VectorMatrix = micm::VectorMatrix<FloatingPointType, 27>;
using Group32VectorMatrix = micm::VectorMatrix<FloatingPointType, 32>;
using Group43VectorMatrix = micm::VectorMatrix<FloatingPointType, 43>;
using Group77VectorMatrix = micm::VectorMatrix<FloatingPointType, 77>;
using Group113VectorMatrix = micm::VectorMatrix<FloatingPointType, 113>;
using Group193VectorMatrix = micm::VectorMatrix<FloatingPointType, 193>;
using Group281VectorMatrix = micm::VectorMatrix<FloatingPointType, 281>;
using Group472VectorMatrix = micm::VectorMatrix<FloatingPointType, 472>;
using Group512VectorMatrix = micm::VectorMatrix<FloatingPointType, 512>;
using Group739VectorMatrix = micm::VectorMatrix<FloatingPointType, 739>;
using Group1130VectorMatrix = micm::VectorMatrix<FloatingPointType, 1130>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group27SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<27>>;
using Group32SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<32>>;
using Group43SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<43>>;
using Group77SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<77>>;
using Group113SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<113>>;
using Group193SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<193>>;
using Group281SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<281>>;
using Group472SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<472>>;
using Group512SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<512>>;
using Group739SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<739>>;
using Group1130SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1130>>;

TEST(RandomCudaProcessSet, Forcing)
{
  testRandomSystemAddForcingTerms<Group1VectorMatrix, Group1CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group3VectorMatrix, Group3CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group27VectorMatrix, Group27CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group32VectorMatrix, Group32CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group43VectorMatrix, Group43CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group77VectorMatrix, Group77CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group113VectorMatrix, Group113CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group193VectorMatrix, Group193CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group281VectorMatrix, Group281CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group472VectorMatrix, Group472CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group512VectorMatrix, Group512CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group739VectorMatrix, Group739CudaDenseMatrix>(404, 500, 400);
  testRandomSystemAddForcingTerms<Group1130VectorMatrix, Group1130CudaDenseMatrix>(404, 500, 400);
}
TEST(RandomCudaProcessSet, Jacobian)
{
  testRandomSystemSubtractJacobianTerms<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      Group1CudaDenseMatrix,
      Group1CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      Group3CudaDenseMatrix,
      Group3CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group27VectorMatrix,
      Group27SparseVectorMatrix,
      Group27CudaDenseMatrix,
      Group27CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group32VectorMatrix,
      Group32SparseVectorMatrix,
      Group32CudaDenseMatrix,
      Group32CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group43VectorMatrix,
      Group43SparseVectorMatrix,
      Group43CudaDenseMatrix,
      Group43CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group77VectorMatrix,
      Group77SparseVectorMatrix,
      Group77CudaDenseMatrix,
      Group77CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group113VectorMatrix,
      Group113SparseVectorMatrix,
      Group113CudaDenseMatrix,
      Group113CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group193VectorMatrix,
      Group193SparseVectorMatrix,
      Group193CudaDenseMatrix,
      Group193CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group281VectorMatrix,
      Group281SparseVectorMatrix,
      Group281CudaDenseMatrix,
      Group281CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group472VectorMatrix,
      Group472SparseVectorMatrix,
      Group472CudaDenseMatrix,
      Group472CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group512VectorMatrix,
      Group512SparseVectorMatrix,
      Group512CudaDenseMatrix,
      Group512CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group739VectorMatrix,
      Group739SparseVectorMatrix,
      Group739CudaDenseMatrix,
      Group739CudaSparseMatrix>(404, 500, 400);
  testRandomSystemSubtractJacobianTerms<
      Group1130VectorMatrix,
      Group1130SparseVectorMatrix,
      Group1130CudaDenseMatrix,
      Group1130CudaSparseMatrix>(404, 500, 400);
}
