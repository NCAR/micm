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
      micm::kokkos::Initialize();
      UpdateDevStruct();
    }

    void UpdateDevStruct()
    {
      micm::kokkos::Initialize();
      devstruct_.number_of_reactants_ = Kokkos::View<std::size_t*>("number_of_reactants", this->number_of_reactants_.size());
      auto h_number_of_reactants = Kokkos::create_mirror_view(devstruct_.number_of_reactants_);
      for (std::size_t i = 0; i < this->number_of_reactants_.size(); ++i)
      {
        h_number_of_reactants(i) = this->number_of_reactants_[i];
      }
      Kokkos::deep_copy(devstruct_.number_of_reactants_, h_number_of_reactants);

      devstruct_.reactant_ids_ = Kokkos::View<std::size_t*>("reactant_ids", this->reactant_ids_.size());
      auto h_reactant_ids = Kokkos::create_mirror_view(devstruct_.reactant_ids_);
      for (std::size_t i = 0; i < this->reactant_ids_.size(); ++i)
      {
        h_reactant_ids(i) = this->reactant_ids_[i];
      }
      Kokkos::deep_copy(devstruct_.reactant_ids_, h_reactant_ids);

      devstruct_.number_of_products_ = Kokkos::View<std::size_t*>("number_of_products", this->number_of_products_.size());
      auto h_number_of_products = Kokkos::create_mirror_view(devstruct_.number_of_products_);
      for (std::size_t i = 0; i < this->number_of_products_.size(); ++i)
      {
        h_number_of_products(i) = this->number_of_products_[i];
      }
      Kokkos::deep_copy(devstruct_.number_of_products_, h_number_of_products);

      devstruct_.product_ids_ = Kokkos::View<std::size_t*>("product_ids", this->product_ids_.size());
      auto h_product_ids = Kokkos::create_mirror_view(devstruct_.product_ids_);
      for (std::size_t i = 0; i < this->product_ids_.size(); ++i)
      {
        h_product_ids(i) = this->product_ids_[i];
      }
      Kokkos::deep_copy(devstruct_.product_ids_, h_product_ids);

      devstruct_.yields_ = Kokkos::View<double*>("yields", this->yields_.size());
      auto h_yields = Kokkos::create_mirror_view(devstruct_.yields_);
      for (std::size_t i = 0; i < this->yields_.size(); ++i)
      {
        h_yields(i) = this->yields_[i];
      }
      Kokkos::deep_copy(devstruct_.yields_, h_yields);

      devstruct_.is_algebraic_variable_ =
          Kokkos::View<uint8_t*>("is_algebraic_variable", this->is_algebraic_variable_.size());
      auto h_is_algebraic_variable = Kokkos::create_mirror_view(devstruct_.is_algebraic_variable_);
      for (std::size_t i = 0; i < this->is_algebraic_variable_.size(); ++i)
      {
        h_is_algebraic_variable(i) = this->is_algebraic_variable_[i];
      }
      Kokkos::deep_copy(devstruct_.is_algebraic_variable_, h_is_algebraic_variable);
    }

    void SetAlgebraicVariableIds(const std::set<std::size_t>& variable_ids)
    {
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetAlgebraicVariableIds(variable_ids);

      devstruct_.is_algebraic_variable_ =
          Kokkos::View<uint8_t*>("is_algebraic_variable", this->is_algebraic_variable_.size());
      auto h_is_algebraic_variable = Kokkos::create_mirror_view(devstruct_.is_algebraic_variable_);
      for (std::size_t i = 0; i < this->is_algebraic_variable_.size(); ++i)
      {
        h_is_algebraic_variable(i) = this->is_algebraic_variable_[i];
      }
      Kokkos::deep_copy(devstruct_.is_algebraic_variable_, h_is_algebraic_variable);
    }

    void SetJacobianFlatIds(const SparseMatrixPolicy& matrix)
    {
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetJacobianFlatIds(matrix);
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

      devstruct_.jacobian_reactant_ids_ =
          Kokkos::View<std::size_t*>("jacobian_reactant_ids", this->jacobian_reactant_ids_.size());
      auto h_jacobian_reactant_ids = Kokkos::create_mirror_view(devstruct_.jacobian_reactant_ids_);
      for (std::size_t i = 0; i < this->jacobian_reactant_ids_.size(); ++i)
      {
        h_jacobian_reactant_ids(i) = this->jacobian_reactant_ids_[i];
      }
      Kokkos::deep_copy(devstruct_.jacobian_reactant_ids_, h_jacobian_reactant_ids);

      devstruct_.jacobian_product_ids_ =
          Kokkos::View<std::size_t*>("jacobian_product_ids", this->jacobian_product_ids_.size());
      auto h_jacobian_product_ids = Kokkos::create_mirror_view(devstruct_.jacobian_product_ids_);
      for (std::size_t i = 0; i < this->jacobian_product_ids_.size(); ++i)
      {
        h_jacobian_product_ids(i) = this->jacobian_product_ids_[i];
      }
      Kokkos::deep_copy(devstruct_.jacobian_product_ids_, h_jacobian_product_ids);

      devstruct_.jacobian_yields_ = Kokkos::View<double*>("jacobian_yields", this->jacobian_yields_.size());
      auto h_jacobian_yields = Kokkos::create_mirror_view(devstruct_.jacobian_yields_);
      for (std::size_t i = 0; i < this->jacobian_yields_.size(); ++i)
      {
        h_jacobian_yields(i) = this->jacobian_yields_[i];
      }
      Kokkos::deep_copy(devstruct_.jacobian_yields_, h_jacobian_yields);

      devstruct_.jacobian_flat_ids_ = Kokkos::View<std::size_t*>("jacobian_flat_ids", this->jacobian_flat_ids_.size());
      auto h_jacobian_flat_ids = Kokkos::create_mirror_view(devstruct_.jacobian_flat_ids_);
      for (std::size_t i = 0; i < this->jacobian_flat_ids_.size(); ++i)
      {
        h_jacobian_flat_ids(i) = this->jacobian_flat_ids_[i];
      }
      Kokkos::deep_copy(devstruct_.jacobian_flat_ids_, h_jacobian_flat_ids);
    }

    void AddForcingTerms(const auto& state, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing) const
    {
      auto d_rate_constants = state.rate_constants_.GetView();
      auto d_state_variables = state_variables.GetView();
      auto d_forcing = forcing.GetView();

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
              double rate = d_rate_constants(rate_idx);
              std::size_t n_reactants = devstruct.number_of_reactants_(i_rxn);
              for (std::size_t i_react = 0; i_react < n_reactants; ++i_react)
              {
                std::size_t reactant_id = devstruct.reactant_ids_(reactant_offset + i_react);
                std::size_t var_idx = (group_id * number_of_species + reactant_id) * L + local_tid;
                rate *= d_state_variables(var_idx);
              }
              for (std::size_t i_react = 0; i_react < n_reactants; ++i_react)
              {
                std::size_t row_id = devstruct.reactant_ids_(reactant_offset + i_react);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t forcing_idx = (group_id * number_of_forcing_species + row_id) * L + local_tid;
                  Kokkos::atomic_add(&d_forcing(forcing_idx), -rate);
                }
              }
              std::size_t n_products = devstruct.number_of_products_(i_rxn);
              for (std::size_t i_prod = 0; i_prod < n_products; ++i_prod)
              {
                std::size_t row_id = devstruct.product_ids_(product_offset + i_prod);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t forcing_idx = (group_id * number_of_forcing_species + row_id) * L + local_tid;
                  Kokkos::atomic_add(&d_forcing(forcing_idx), devstruct.yields_(product_offset + i_prod) * rate);
                }
              }
              reactant_offset += n_reactants;
              product_offset += n_products;
            }
          });
    }

    void SubtractJacobianTerms(const auto& state, const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian)
        const
    {
      auto d_rate_constants = state.rate_constants_.GetView();
      auto d_state_variables = state_variables.GetView();
      auto d_jacobian = jacobian.GetView();
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
              double d_rate_d_ind = d_rate_constants(rate_idx);
              for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
              {
                std::size_t species_id = devstruct.jacobian_reactant_ids_(reactant_offset + i_react);
                std::size_t var_idx = (group_id * number_of_species + species_id) * L + local_tid;
                d_rate_d_ind *= d_state_variables(var_idx);
              }
              for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_; ++i_dep)
              {
                std::size_t row_id = devstruct.jacobian_reactant_ids_(reactant_offset + i_dep);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t jac_idx = jac_group_offset + devstruct.jacobian_flat_ids_(flat_id_offset) + local_tid;
                  Kokkos::atomic_add(&d_jacobian(jac_idx), d_rate_d_ind);
                }
                flat_id_offset++;
              }
              if (!devstruct.is_algebraic_variable_(process_info.independent_id_))
              {
                std::size_t jac_idx = jac_group_offset + devstruct.jacobian_flat_ids_(flat_id_offset) + local_tid;
                Kokkos::atomic_add(&d_jacobian(jac_idx), d_rate_d_ind);
              }
              flat_id_offset++;
              for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
              {
                std::size_t row_id = devstruct.jacobian_product_ids_(product_offset + i_dep);
                if (!devstruct.is_algebraic_variable_(row_id))
                {
                  std::size_t jac_idx = jac_group_offset + devstruct.jacobian_flat_ids_(flat_id_offset) + local_tid;
                  Kokkos::atomic_sub(
                      &d_jacobian(jac_idx), devstruct.jacobian_yields_(product_offset + i_dep) * d_rate_d_ind);
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
