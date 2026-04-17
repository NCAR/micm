// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <Kokkos_Core.hpp>
#include <cstddef>
#include <cstdint>

namespace micm
{
  namespace kokkos
  {
    /// This struct holds information about a process for the Jacobian calculation
    struct ProcessInfoParam
    {
      std::size_t process_id_;
      std::size_t independent_id_;
      std::size_t number_of_dependent_reactants_;
      std::size_t number_of_products_;
    };

    /// This struct holds the (1) pointer to, and (2) size of
    ///   each constant data member from the class "ProcessSet";
    /// This struct could be used within Kokkos kernels;
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
