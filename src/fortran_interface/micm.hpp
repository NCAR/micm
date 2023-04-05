#pragma once

#include <micm/solver/solver.hpp>
#include <micm/solver/state.hpp>

namespace micm {

  class MICM {
    Solver* solver_;
    State* state_;
  };
}