// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <Kokkos_Core.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace micm
{
  namespace kokkos
  {
    /// @brief Copy a host std::vector into a new Kokkos::View on the default execution space
    template<typename T>
    Kokkos::View<T*> CopyVectorToView(const std::string& name, const std::vector<T>& src)
    {
      Kokkos::View<T*> view(name, src.size());
      auto host = Kokkos::create_mirror_view(view);
      for (std::size_t i = 0; i < src.size(); ++i)
      {
        host(i) = src[i];
      }
      Kokkos::deep_copy(view, host);
      return view;
    }

    /// This struct holds information about a process for the Jacobian calculation
    struct ProcessInfoParam
    {
      std::size_t process_id_;
      std::size_t independent_id_;
      std::size_t number_of_dependent_reactants_;
      std::size_t number_of_products_;
    };

    /// This struct holds device-side copies of ProcessSet data for use in Kokkos kernels
    struct ProcessSetParam
    {
      Kokkos::View<std::size_t*> number_of_reactants_;
      Kokkos::View<std::size_t*> reactant_ids_;
      Kokkos::View<std::size_t*> number_of_products_;
      Kokkos::View<std::size_t*> product_ids_;
      Kokkos::View<double*> yields_;
      Kokkos::View<ProcessInfoParam*> jacobian_process_info_;
      Kokkos::View<std::size_t*> jacobian_reactant_ids_;
      Kokkos::View<std::size_t*> jacobian_product_ids_;
      Kokkos::View<double*> jacobian_yields_;
      Kokkos::View<std::size_t*> jacobian_flat_ids_;
      Kokkos::View<uint8_t*> is_algebraic_variable_;
    };
  }  // namespace kokkos
}  // namespace micm
