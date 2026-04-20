// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_param.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/process/process_set.hpp>

#include <Kokkos_Core.hpp>

namespace micm
{
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  class KokkosProcessSet : public ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>
  {
   public:
    kokkos::ProcessSetParam devstruct_;

    KokkosProcessSet() = default;

    KokkosProcessSet(const std::vector<Process>& processes, const std::unordered_map<std::string, std::size_t>& variable_map)
        : ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>(processes, variable_map)
    {
      UpdateDevStruct();
    }

    void UpdateDevStruct()
    {
      devstruct_.number_of_reactants_ = kokkos::CopyVectorToView("number_of_reactants", this->number_of_reactants_);
      devstruct_.reactant_ids_ = kokkos::CopyVectorToView("reactant_ids", this->reactant_ids_);
      devstruct_.number_of_products_ = kokkos::CopyVectorToView("number_of_products", this->number_of_products_);
      devstruct_.product_ids_ = kokkos::CopyVectorToView("product_ids", this->product_ids_);
      devstruct_.yields_ = kokkos::CopyVectorToView("yields", this->yields_);

      // Convert bool vector to uint8_t for device compatibility
      std::vector<uint8_t> algebraic_flags(this->is_algebraic_variable_.begin(), this->is_algebraic_variable_.end());
      devstruct_.is_algebraic_variable_ = kokkos::CopyVectorToView("is_algebraic_variable", algebraic_flags);
    }

    void SetAlgebraicVariableIds(const std::set<std::size_t>& variable_ids)
    {
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetAlgebraicVariableIds(variable_ids);

      std::vector<uint8_t> algebraic_flags(this->is_algebraic_variable_.begin(), this->is_algebraic_variable_.end());
      devstruct_.is_algebraic_variable_ = kokkos::CopyVectorToView("is_algebraic_variable", algebraic_flags);
    }

    void SetJacobianFlatIds(const SparseMatrixPolicy& matrix)
    {
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetJacobianFlatIds(matrix);

      // Copy Jacobian process info (struct-of-arrays style)
      devstruct_.jacobian_process_info_ =
          Kokkos::View<kokkos::ProcessInfoParam*>("jacobian_process_info", this->jacobian_process_info_.size());
      auto h_jacobian_process_info = Kokkos::create_mirror_view(devstruct_.jacobian_process_info_);
      for (std::size_t i = 0; i < this->jacobian_process_info_.size(); ++i)
      {
        h_jacobian_process_info(i).process_id_ = this->jacobian_process_info_[i].process_id_;
        h_jacobian_process_info(i).independent_id_ = this->jacobian_process_info_[i].independent_id_;
        h_jacobian_process_info(i).number_of_dependent_reactants_ =
            this->jacobian_process_info_[i].number_of_dependent_reactants_;
        h_jacobian_process_info(i).number_of_products_ = this->jacobian_process_info_[i].number_of_products_;
      }
      Kokkos::deep_copy(devstruct_.jacobian_process_info_, h_jacobian_process_info);

      devstruct_.jacobian_reactant_ids_ = kokkos::CopyVectorToView("jacobian_reactant_ids", this->jacobian_reactant_ids_);
      devstruct_.jacobian_product_ids_ = kokkos::CopyVectorToView("jacobian_product_ids", this->jacobian_product_ids_);
      devstruct_.jacobian_yields_ = kokkos::CopyVectorToView("jacobian_yields", this->jacobian_yields_);
      devstruct_.jacobian_flat_ids_ = kokkos::CopyVectorToView("jacobian_flat_ids", this->jacobian_flat_ids_);
    }

    template<typename StatePolicy>
    void AddForcingTerms(const StatePolicy& state, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing)
        const
    {
      auto view_rate_constants = state.rate_constants_.GetView();
      auto view_state_variables = state_variables.GetView();
      auto view_forcing = forcing.GetView();
      auto devstruct = this->devstruct_;

      std::size_t number_of_grid_cells = state_variables.NumRows();
      std::size_t number_of_reactions = this->number_of_reactants_.size();
      std::size_t number_of_species = state_variables.NumColumns();
      std::size_t number_of_forcing_species = forcing.NumColumns();
      std::size_t L = DenseMatrixPolicy::GroupVectorSize();

      Kokkos::parallel_for(
          "AddForcingTerms", number_of_grid_cells, KOKKOS_LAMBDA(const std::size_t i_cell) {
            std::size_t group_id = i_cell / L;
            std::size_t local_tid = i_cell % L;

            std::size_t reactant_offset = 0;
            std::size_t product_offset = 0;
            for (std::size_t i_rxn = 0; i_rxn < number_of_reactions; ++i_rxn)
            {
              std::size_t rate_idx = (group_id * number_of_reactions + i_rxn) * L + local_tid;
              double rate = view_rate_constants(rate_idx);
              std::size_t n_reactants = devstruct.number_of_reactants_(i_rxn);
              for (std::size_t i_react = 0; i_react < n_reactants; ++i_react)
              {
                std::size_t reactant_id = devstruct.reactant_ids_(reactant_offset + i_react);
                std::size_t var_idx = (group_id * number_of_species + reactant_id) * L + local_tid;
                rate *= view_state_variables(var_idx);
              }
              for (std::size_t i_react = 0; i_react < n_reactants; ++i_react)
              {
                std::size_t row_id = devstruct.reactant_ids_(reactant_offset + i_react);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t forcing_idx = (group_id * number_of_forcing_species + row_id) * L + local_tid;
                  view_forcing(forcing_idx) -= rate;
                }
              }
              std::size_t n_products = devstruct.number_of_products_(i_rxn);
              for (std::size_t i_prod = 0; i_prod < n_products; ++i_prod)
              {
                std::size_t row_id = devstruct.product_ids_(product_offset + i_prod);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t forcing_idx = (group_id * number_of_forcing_species + row_id) * L + local_tid;
                  view_forcing(forcing_idx) += devstruct.yields_(product_offset + i_prod) * rate;
                }
              }
              reactant_offset += n_reactants;
              product_offset += n_products;
            }
          });
    }

    template<typename StatePolicy>
    void SubtractJacobianTerms(
        const StatePolicy& state,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const
    {
      auto view_rate_constants = state.rate_constants_.GetView();
      auto view_state_variables = state_variables.GetView();
      auto view_jacobian = jacobian.GetView();
      auto devstruct = this->devstruct_;

      std::size_t number_of_grid_cells = state_variables.NumRows();
      std::size_t number_of_species = state_variables.NumColumns();
      std::size_t number_of_process_infos = devstruct.jacobian_process_info_.extent(0);
      std::size_t number_of_reactions = this->number_of_reactants_.size();
      std::size_t number_of_non_zeros = jacobian.FlatBlockSize();
      std::size_t L = DenseMatrixPolicy::GroupVectorSize();
      std::size_t jacobian_group_size = number_of_non_zeros * L;

      Kokkos::parallel_for(
          "SubtractJacobianTerms", number_of_grid_cells, KOKKOS_LAMBDA(const std::size_t i_cell) {
            std::size_t group_id = i_cell / L;
            std::size_t local_tid = i_cell % L;
            std::size_t jac_group_offset = group_id * jacobian_group_size;

            std::size_t reactant_offset = 0;
            std::size_t product_offset = 0;
            std::size_t flat_id_offset = 0;
            for (std::size_t i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
            {
              const auto& process_info = devstruct.jacobian_process_info_(i_proc);
              std::size_t rate_idx = (group_id * number_of_reactions + process_info.process_id_) * L + local_tid;
              double d_rate_d_ind = view_rate_constants(rate_idx);
              for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
              {
                std::size_t species_id = devstruct.jacobian_reactant_ids_(reactant_offset + i_react);
                std::size_t var_idx = (group_id * number_of_species + species_id) * L + local_tid;
                d_rate_d_ind *= view_state_variables(var_idx);
              }
              for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_; ++i_dep)
              {
                std::size_t row_id = devstruct.jacobian_reactant_ids_(reactant_offset + i_dep);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t jac_idx = jac_group_offset + devstruct.jacobian_flat_ids_(flat_id_offset) + local_tid;
                  view_jacobian(jac_idx) += d_rate_d_ind;
                }
                flat_id_offset++;
              }
              if (!devstruct.is_algebraic_variable_(process_info.independent_id_))
              {
                std::size_t jac_idx = jac_group_offset + devstruct.jacobian_flat_ids_(flat_id_offset) + local_tid;
                view_jacobian(jac_idx) += d_rate_d_ind;
              }
              flat_id_offset++;
              for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
              {
                std::size_t row_id = devstruct.jacobian_product_ids_(product_offset + i_dep);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t jac_idx = jac_group_offset + devstruct.jacobian_flat_ids_(flat_id_offset) + local_tid;
                  view_jacobian(jac_idx) -= devstruct.jacobian_yields_(product_offset + i_dep) * d_rate_d_ind;
                }
                flat_id_offset++;
              }
              reactant_offset += process_info.number_of_dependent_reactants_;
              product_offset += process_info.number_of_products_;
            }
          });
    }
  };
}  // namespace micm
