#include "../../util/test_sparse_matrix_policy.hpp"
#include "cuda_matrix_utils.cuh"

#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(CudaSparseMatrix, ZeroMatrix)
{
  auto matrix = testZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(CudaSparseMatrix, ConstZeroMatrix)
{
  const auto matrix = testConstZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(CudaSparseMatrix, CopyAssignmentZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  EXPECT_EQ(matrix.AsVector().size(), 2);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, CopyAssignmentConstZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  const micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  auto oneMatrix = matrix;
  oneMatrix.CopyToDevice();

  auto param = oneMatrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  EXPECT_EQ(matrix.AsVector().size(), 2);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, CopyAssignmentDeSynchedHostZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  matrix.CopyToDevice();

  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  auto oneMatrix = matrix;

  EXPECT_EQ(2, matrix.AsVector().size());
  EXPECT_EQ(2, oneMatrix.AsVector().size());

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, MoveAssignmentConstZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  auto oneMatrix = std::move(matrix);
  if (matrix.AsDeviceParam().d_data_ != nullptr)
  {
    throw std::runtime_error(
        "The 'd_data_' pointer of oneMatrix is not initialized to a null pointer in the move constructor.");
  }

  EXPECT_EQ(2, oneMatrix.AsVector().size());
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToDevice();
  auto param = oneMatrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, MoveAssignmentDeSyncedHostZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  EXPECT_EQ(2, matrix.AsVector().size());
  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  auto oneMatrix = std::move(matrix);
  if (matrix.AsDeviceParam().d_data_ != nullptr)
  {
    throw std::runtime_error(
        "The 'd_data_' pointer of oneMatrix is not initialized to a null pointer in the move constructor.");
  }

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, SetScalar)
{
  testSetScalar<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(CudaSparseMatrix, SingleBlockMatrix)
{
  auto matrix = testSingleBlockMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<4>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 16);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[16], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 12);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[12], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(matrix.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(CudaSparseMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<2>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 8);
    EXPECT_EQ(matrix.AsVector()[8], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(CudaSparseMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<2>>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 6);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 2, 1);
    EXPECT_EQ(elem, 14);
    matrix.AsVector()[elem] = 31;
    EXPECT_EQ(matrix.AsVector()[14], 31);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(4), 2);
}
