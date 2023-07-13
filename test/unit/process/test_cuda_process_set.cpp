#include <gtest/gtest.h>
#include <micm/process/process_set_cuda.cuh>
#include <micm/process/process_set.hpp>
#include <micm/util/vector_matrix.hpp>
#include <iostream>
#include <random>
#include <chrono>
#include <functional>
#include <vector>

using yields = std::pair<micm::Species, double>;
using index_pair = std::pair<std::size_t, std::size_t>;

void compare_pair(const index_pair& a, const index_pair& b)
{
  EXPECT_EQ(a.first, b.first);
  EXPECT_EQ(a.second, b.second);
}



template<template<class> class MatrixPolicy>
void testRandomSystem(std::size_t n_cells, std::size_t n_reactions, std::size_t n_species)
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
  micm::State state{ micm::StateParameters{ .state_variable_names_{ species_names },
                                                          .number_of_grid_cells_ = n_cells,
                                                          .number_of_custom_parameters_ = 0,
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

  micm::ProcessSet set{ processes, state };

  for (auto& elem : state.variables_.AsVector())
    elem = get_double();

  micm::Matrix <double> rate_constants{ n_cells, n_reactions };
  for (auto& elem : rate_constants.AsVector())
    elem = get_double();

  micm::Matrix <double> cpu_forcing{ n_cells, n_species, 1000.0};
  micm::Matrix <double> gpu_forcing{ };
  gpu_forcing = cpu_forcing; 

  
  //kernel function call 
  // double t0 = 0.0; 
  // for (int i = 0; i < 100; i++)
  // {
  //   auto start = std::chrono::steady_clock::now(); 
  //   set.CudaAddForcingTerms(rate_constants, state.variables_, gpu_forcing); 
  //   auto end = std::chrono::steady_clock::now(); 
  //   std::chrono::duration<double> duration = end - start;
  //   t0 = t0 + duration.count(); 
  // }
  // std::cout << "time performance: "<< t0/100 <<std::endl; 

  //kernel function call 
  set.CudaAddForcingTerms(rate_constants, state.variables_, gpu_forcing); 
    
  //CPU function call
  set.AddForcingTerms(rate_constants, state.variables_, cpu_forcing); 

  //checking accuracy with comparison between CPU and GPU result 
  std::vector<double>cpu_forcing_vector = cpu_forcing.AsVector(); 
  std::vector<double>gpu_forcing_vector = gpu_forcing.AsVector(); 

  for (int i = 0; i < cpu_forcing_vector.size(); i++){
    double a = cpu_forcing_vector[i];
    double b = gpu_forcing_vector[i];
    EXPECT_NEAR(a, b, std::abs(a+b)*1.0e-9);
 }
}

template<class T>
using Group1000VectorMatrix = micm::VectorMatrix<T, 1000>;
template<class T>
using Group10000VectorMatrix = micm::VectorMatrix<T, 10000>;
template<class T>
using Group100000VectorMatrix = micm::VectorMatrix<T, 100000>;
template<class T>
using Group1000000VectorMatrix = micm::VectorMatrix<T, 1000000>;

TEST(RandomProcessSet, Matrix)
{
  std::cout << "system with 500 reactions and 400 species"<<std::endl; 
  testRandomSystem<Group1000VectorMatrix>(1000, 500, 400);
  testRandomSystem<Group10000VectorMatrix>(10000, 500, 400);
  testRandomSystem<Group100000VectorMatrix>(100000, 500, 400);
  testRandomSystem<Group1000000VectorMatrix>(1000000, 500, 400);

  std::cout << "system with 100 reactions and 80 species"<<std::endl; 
  testRandomSystem<Group1000VectorMatrix>(1000, 100, 80);
  testRandomSystem<Group10000VectorMatrix>(10000, 100, 80);
  testRandomSystem<Group100000VectorMatrix>(100000, 100, 80);
  testRandomSystem<Group1000000VectorMatrix>(1000000, 100, 80);
}
