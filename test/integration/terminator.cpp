#include "terminator.hpp"

#include <gtest/gtest.h>

#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

TEST(RosenbrockSolver, Terminator)
{
  TestTerminator<micm::Matrix, micm::SparseMatrix, micm::LinearSolver<double, micm::SparseMatrix>>(1);
  TestTerminator<micm::Matrix, micm::SparseMatrix, micm::LinearSolver<double, micm::SparseMatrix>>(2);
  TestTerminator<micm::Matrix, micm::SparseMatrix, micm::LinearSolver<double, micm::SparseMatrix>>(3);
  TestTerminator<micm::Matrix, micm::SparseMatrix, micm::LinearSolver<double, micm::SparseMatrix>>(4);
}

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

TEST(RosenbrockSolver, VectorTerminator)
{
  TestTerminator<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(1);
  TestTerminator<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(4);
  TestTerminator<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(3);
  TestTerminator<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(2);
}