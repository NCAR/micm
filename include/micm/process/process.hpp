/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

namespace micm
{

  class NoPolicy {
   private:
   public:
  };

  class GenerateForcingPolicy {
   private:
   public:
  };

  class GenerateJacobianPolicy {
   private:
   public:
  };

  template<
    class ReactionType = NoPolicy,
    class ForcingPolicy = NoPolicy,
    class JacobianPolicy = NoPolicy
  >
  class Process
  {
    ReactionType reaction_type_;
    ForcingPolicy forcing_policy_;
    JacobianPolicy jacobian_policy_;
   private:
   public:

  };

}  // namespace micm
