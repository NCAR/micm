#include <micm/system/species.hpp>

namespace micm
{

  Species::Species() = default;

  Species::Species(std::string name)
      : name_(std::move(name))
  {
  }

  Species::Species(std::string name, std::vector<Property<double>> properties)
      : name_(std::move(name)),
        properties_(std::move(properties))
  {
  }

}  // namespace micm