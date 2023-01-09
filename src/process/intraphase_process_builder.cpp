#include <micm/process/intraphase_process_builder.hpp>

namespace micm
{

  IntraPhaseProcessBuilder& micm::IntraPhaseProcessBuilder::For(const Phase& phase)
  {
    return *this;
  }
  IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::Reacting(const Species& reactant)
  {
    return *this;
    // TODO: insert return statement here
  }

  IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::Producing(const Species& product)
  {
    return *this;
  }

  IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::WithYield(const double yield)
  {
    return *this;
  }

  IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::WithRateConstant(const RateConstant& rate_constant)
  {
    return *this;
  }

  IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::Build()
  {
    return *this;
  }

}  // namespace micm