#include <micm/system/condition.hpp>

namespace micm
{

  Condition::Condition() = default;
  Condition::Condition(std::string name, std::string units)
      : name_(std::move(name)),
        units_(std::move(units)){};

}  // namespace micm