#include <micm/system/condition.hpp>

namespace micm
{

  Condition::Condition() = default;
  Condition::Condition(const std::string& name, const std::string& units) : name_(name), units_(units){};

}  // namespace micm