#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <vector>

namespace micm
{
  void AddForcingTerms_kernelSetup(
      const Matrix<double>& rate_constants,
      const Matrix<double>& state_variables,
      Matrix<double>& forcing,
      const std::vector<std::size_t>& number_of_reactants_,
      const std::vector<std::size_t>& reactant_ids_,
      const std::vector<std::size_t>& number_of_products_,
      const std::vector<std::size_t>& product_ids_,
      const std::vector<double>& yields_,
      const std::vector<std::size_t>& jacobian_flat_ids_);
}  // namespace micm
