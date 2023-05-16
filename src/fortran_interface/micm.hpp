#pragma once

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>

namespace micm {

  class MICM {
    public:
      RosenbrockSolver* solver_;
      State* state_;
  };

}