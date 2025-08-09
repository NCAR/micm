#pragma once

#include <memory>
#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>  // Assuming this is where the base is
#include <micm/system/conditions.hpp>

namespace micm
{
  class TestTransferCoefficient : public TransferCoefficient
  {
  public:
    TestTransferCoefficient() = default;

    virtual ~TestTransferCoefficient() = default;

    /// @brief Clone this object (returns a new copy on the heap)
    virtual std::unique_ptr<TransferCoefficient> Clone() const override
    {
      return std::make_unique<TestTransferCoefficient>(*this);
    }

    /// @brief Simple constant value for testing
    virtual double Calculate() const override
    {
      return 1.0;
    }

    /// @brief Simple override for environmental condition-based calculation
    virtual double Calculate(const Conditions& conditions) const override
    {
      // Just return a dummy value for testing
      return 1.0;
    }
  };
}
