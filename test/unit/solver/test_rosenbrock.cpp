#include <gtest/gtest.h>

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

TEST(ChapmanODESolver, DefaultConstructor)
{
  micm::RosenbrockSolver<micm::Matrix, micm::SparseMatrix> solver{};
}