#include <micm/process/process_set.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <algorithm>

namespace
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");

}  // namespace
