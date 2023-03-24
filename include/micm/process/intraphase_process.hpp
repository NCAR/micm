/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/species.hpp>
#include <vector>

namespace micm
{

  /**
   * @brief An intraphase process
   *
   */
  template<class Rate>
  class IntraPhaseProcess
  {
   public:
    std::vector<Species> reactants_;
    std::vector<Species> products_;
    Rate rate_;

   public:
    /// @brief Default constructor
    IntraPhaseProcess();

    /// @brief
    /// @param reactants
    /// @param products
    /// @param rate
    IntraPhaseProcess(std::vector<Species> reactants, std::vector<Species> products, Rate rate);

    /// @brief Copy constructor
    IntraPhaseProcess(const IntraPhaseProcess& other);

    /// @brief Move constructor
    IntraPhaseProcess(IntraPhaseProcess&& other);

    /// @brief Copy assignment operator
    IntraPhaseProcess& operator=(IntraPhaseProcess& other);

    /// @brief Move assignment operator
    IntraPhaseProcess& operator=(IntraPhaseProcess&& other);
  };

  template<class Rate>
  inline IntraPhaseProcess<Rate>::IntraPhaseProcess()
      : reactants_(),
        products_(),
        rate_()
  {
  }

  template<class Rate>
  inline IntraPhaseProcess<Rate>::IntraPhaseProcess(std::vector<Species> reactants, std::vector<Species> products, Rate rate)
      : reactants_(reactants),
        products_(products),
        rate_(rate)
  {
  }

  template<class Rate>
  inline IntraPhaseProcess<Rate>::IntraPhaseProcess(const IntraPhaseProcess<Rate>& other)
      : reactants_(other.reactants_),
        products_(other.products_),
        rate_(other.rate_)
  {
  }

  template<class Rate>
  inline IntraPhaseProcess<Rate>::IntraPhaseProcess(IntraPhaseProcess<Rate>&& other)
      : reactants_(std::move(other.reactants_)),
        products_(std::move(other.products_)),
        rate_(std::move(other.rate_))
  {
  }

  template<class Rate>
  inline IntraPhaseProcess<Rate>& IntraPhaseProcess<Rate>::operator=(IntraPhaseProcess<Rate>& other)
  {
    if (this != &other)
    {
      reactants_ = other.reactants_;
      products_ = other.products_;
      rate_ = other.rate_;
    }
    return *this;
  }

  template<class Rate>
  inline IntraPhaseProcess<Rate>& IntraPhaseProcess<Rate>::operator=(IntraPhaseProcess<Rate>&& other)
  {
    if (this != &other)
    {
      reactants_ = std::move(other.reactants_);
      products_ = std::move(other.products_);
      rate_ = std::move(other.rate_);
    }
    return *this;
  }

}  // namespace micm