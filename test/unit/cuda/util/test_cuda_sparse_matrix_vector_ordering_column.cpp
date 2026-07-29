#include "../../util/test_sparse_matrix_policy.hpp"
#include "cuda_matrix_utils.cuh"

#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering_compressed_sparse_column.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>

/* Below are the policy Tests on the CPU */

TEST(SparseVectorCompressedColumnMatrix, ZeroMatrix)
{
  TestZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstZeroMatrix)
{
  TestConstZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, SetScalar)
{
  TestSetScalar<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, AddToDiagonal)
{
  TestAddToDiagonal<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestAddToDiagonal<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestAddToDiagonal<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, Print)
{
  TestPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, PrintNonZero)
{
  TestPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

/* These are the customized Tests for GPU only */

TEST(CudaSparseMatrix, CopyAssignmentZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

  matrix.CopyToDevice();

  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  auto oneMatrix = matrix;

  EXPECT_EQ(2, matrix.AsVector().size());
  EXPECT_EQ(2, oneMatrix.AsVector().size());

  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToHost();

  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, CopyAssignmentConstZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  const micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

  auto oneMatrix = matrix;
  oneMatrix.CopyToDevice();

  auto param = oneMatrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  EXPECT_EQ(matrix.AsVector().size(), 2);

  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, MoveAssignmentConstZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

  auto oneMatrix = std::move(matrix);
  if (matrix.AsDeviceParam().d_data_ != nullptr)  // NOLINT(bugprone-use-after-move): checks moved-from state
  {
    throw std::runtime_error(
        "The 'd_data_' pointer of oneMatrix is not initialized to a null pointer in the move constructor.");
  }

  EXPECT_EQ(2, oneMatrix.AsVector().size());
  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToDevice();
  auto param = oneMatrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseMatrix, MoveAssignmentZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<micm::Real, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

  EXPECT_EQ(2, matrix.AsVector().size());
  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (const auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  auto oneMatrix = std::move(matrix);
  if (matrix.AsDeviceParam().d_data_ != nullptr)  // NOLINT(bugprone-use-after-move): checks moved-from state
  {
    throw std::runtime_error(
        "The 'd_data_' pointer of oneMatrix is not initialized to a null pointer in the move constructor.");
  }

  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (const auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

template<micm::Index cuda_matrix_vector_length>
void TestSingleBlockMatrixAddOneElement()
{
  auto matrix = TestSingleBlockMatrix<
      micm::CudaSparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<cuda_matrix_vector_length>>();

  {
    micm::Index elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 3 * cuda_matrix_vector_length);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[elem], 42);
  }
  {
    micm::Index elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 4 * cuda_matrix_vector_length);
    matrix.AsVector()[elem] = 39;
    EXPECT_EQ(matrix.AsVector()[elem], 39);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), cuda_matrix_vector_length);
  EXPECT_EQ(matrix.GroupSize(), 5 * cuda_matrix_vector_length);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  micm::Index elem_id = 2;  // in this example, 2 refers to matrix[2][1] which is non-zero
  micm::cuda::SparseMatrixAddOneElementDriver(param, elem_id, 0, cuda_matrix_vector_length);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][2][1], 46);
}

TEST(CudaSparseMatrix, SingleBlockMatrixAddOneElement)
{
  TestSingleBlockMatrixAddOneElement<3>();
  TestSingleBlockMatrixAddOneElement<32>();
  TestSingleBlockMatrixAddOneElement<37>();
  TestSingleBlockMatrixAddOneElement<65>();
}

template<micm::Index cuda_matrix_vector_length>
void TestMultiBlockMatrixAddOneElement()
{
  auto matrix = TestMultiBlockMatrix<
      micm::CudaSparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<cuda_matrix_vector_length>>();

  {
    micm::Index elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 4 * cuda_matrix_vector_length);
    auto idx = 40;
    elem = matrix.VectorIndex(idx, 2, 3);
    auto a = std::floor(static_cast<micm::Real>(idx) / cuda_matrix_vector_length);
    auto b = matrix.VectorIndex(0, 2, 3) / cuda_matrix_vector_length;
    auto c = idx % cuda_matrix_vector_length;
    EXPECT_EQ(elem, 5 * cuda_matrix_vector_length * a + b * cuda_matrix_vector_length + c);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), cuda_matrix_vector_length);
  EXPECT_EQ(matrix.GroupSize(), cuda_matrix_vector_length * 5);
  // Work around NVHPC compiler bug - force runtime evaluation
  micm::Real num_cells_d = 53.0;
  auto vec_length_d = static_cast<micm::Real>(cuda_matrix_vector_length);
  auto expected_groups = static_cast<micm::Index>(std::ceil(num_cells_d / vec_length_d));
  EXPECT_EQ(matrix.NumberOfGroups(53), expected_groups);

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  micm::Index elem_id = 4;  // in this example, 4 refers to matrix[2,3] which is non-zero
  micm::Index grid_id = 33;
  micm::cuda::SparseMatrixAddOneElementDriver(param, elem_id, grid_id, cuda_matrix_vector_length);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[grid_id][2][3], 80);
}

TEST(CudaSparseMatrix, MultiBlockMatrixAddOneElement)
{
  TestMultiBlockMatrixAddOneElement<3>();
  TestMultiBlockMatrixAddOneElement<32>();
  TestMultiBlockMatrixAddOneElement<37>();
  TestMultiBlockMatrixAddOneElement<65>();
}
