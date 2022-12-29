#include <micm/system/species.hpp>

namespace micm
{

  Species::Species() = default;

  Species::Species(std::string name) : name_(name)
  {
  }

  Species::Species(std::string name, std::vector<Property<double>> properties) : name_(name), properties_(properties)
  {
  }

}  // namespace micm