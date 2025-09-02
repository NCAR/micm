// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process_error.hpp>
#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <functional>
#include <iostream>
#include <vector>

// This is a design concept to set forcing terms, and later Jacobian terms,
// for phase transfer processes that involve aerosol models.
// This structure enables a callback-based interaction between two classes:
// - PhaseTransferProcessSet
// - external classes such as AerosolModel
// PhaseTransferProcessSet owns the logic to assemble and manage forcing terms.
// It allows external classes to register callbacks via RegisterForcingCallback().
// When AddForcingTerms(..) is called, the callback is triggered to allow the model
// to modify the forcing terms.
// An external model has full control over how it modifies the forcing data, and
// PhaseTransferProcessSet doesn't need to know any of its internal logic.
namespace micm
{
  class PhaseTransferProcessSet
  {
   public:
    using CallbackType = std::function<void(std::vector<std::vector<double>>& forcing)>;

    void RegisterForcingCallback(CallbackType callback_type)
    {
      callback_ = callback_type;
    }

    void AddForcingTerms(std::vector<std::vector<double>>& forcing)
    {
      // TODO (jiwon) - Set the model's constraint contributions for forcing
      if (callback_)
      {
        std::cout << "[PhaseTransferProcessSet] Calling a model to set forcing terms" << std::endl;
        callback_(forcing);
      }

      // TODO (jiwon) - Set the forcing terms
      constexpr double TEMP_VALUE_FOR_TEST = 0.3333;
      forcing[0][0] = TEMP_VALUE_FOR_TEST;
    }

   private:
    CallbackType callback_;
  };

}  // namespace micm
