#include <micm/system/phase.hpp>

namespace micm
{

  Phase::Phase() = default;

  Phase::Phase(const Phase& other)
  {
  }

  Phase::Phase(std::vector<Species> species)
    : species_(std::move(species))
  {
  }

  Phase& Phase::operator=(const Phase& other)
  {
    // TODO: insert return statement here
  }

  size_t Phase::Size()
  {
    // TODO
  }
}  // namespace micm