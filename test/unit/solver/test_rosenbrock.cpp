#include <gtest/gtest.h>

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

TEST(ChapmanODESolver, DefaultConstructor)
{
  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{};
}