#include <micm/process/intraphase_process_builder.hpp>

namespace micm
{

  IntraPhaseProcessBuilder& micm::IntraPhaseProcessBuilder::For(const Phase& phase)
  {
    return *this;
  }
  IntraPhaseProcessBuilder& micm::IntraPhaseProcessBuilder::With(const Species& phase)
  {
    return *this;
  }

  IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::Build()
  {
    return *this;
  }

}  // namespace micm