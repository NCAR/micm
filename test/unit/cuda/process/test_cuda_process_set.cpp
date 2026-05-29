#include "../../process/test_process_set_policy.hpp"

#include <micm/cuda/process/cuda_process_set.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

using FloatingPointType = double;

using Group1CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 1>;
using Group3CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 3>;
using Group27CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 27>;
using Group32CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 32>;
using Group43CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 43>;
using Group77CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 77>;
using Group113CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 113>;
using Group193CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 193>;
using Group281CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 281>;
using Group472CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 472>;
using Group512CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 512>;
using Group739CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 739>;
using Group1130CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 1130>;

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group3CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group27CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<27>>;
using Group32CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<32>>;
using Group43CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<43>>;
using Group77CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<77>>;
using Group113CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<113>>;
using Group193CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<193>>;
using Group281CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<281>>;
using Group472CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<472>>;
using Group512CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<512>>;
using Group739CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<739>>;
using Group1130CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1130>>;

/* These are the policy tests on the GPU */

TEST(RandomCudaProcessSet, CudaMatrix)
{
  TestRandomSystem<
      Group1CudaDenseMatrix,
      Group1CudaSparseMatrix,
      micm::CudaProcessSet<Group1CudaDenseMatrix, Group1CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group3CudaDenseMatrix,
      Group3CudaSparseMatrix,
      micm::CudaProcessSet<Group3CudaDenseMatrix, Group3CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group27CudaDenseMatrix,
      Group27CudaSparseMatrix,
      micm::CudaProcessSet<Group27CudaDenseMatrix, Group27CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group32CudaDenseMatrix,
      Group32CudaSparseMatrix,
      micm::CudaProcessSet<Group32CudaDenseMatrix, Group32CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group43CudaDenseMatrix,
      Group43CudaSparseMatrix,
      micm::CudaProcessSet<Group43CudaDenseMatrix, Group43CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group77CudaDenseMatrix,
      Group77CudaSparseMatrix,
      micm::CudaProcessSet<Group77CudaDenseMatrix, Group77CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group113CudaDenseMatrix,
      Group113CudaSparseMatrix,
      micm::CudaProcessSet<Group113CudaDenseMatrix, Group113CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group193CudaDenseMatrix,
      Group193CudaSparseMatrix,
      micm::CudaProcessSet<Group193CudaDenseMatrix, Group193CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group281CudaDenseMatrix,
      Group281CudaSparseMatrix,
      micm::CudaProcessSet<Group281CudaDenseMatrix, Group281CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group472CudaDenseMatrix,
      Group472CudaSparseMatrix,
      micm::CudaProcessSet<Group472CudaDenseMatrix, Group472CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group512CudaDenseMatrix,
      Group512CudaSparseMatrix,
      micm::CudaProcessSet<Group512CudaDenseMatrix, Group512CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group739CudaDenseMatrix,
      Group739CudaSparseMatrix,
      micm::CudaProcessSet<Group739CudaDenseMatrix, Group739CudaSparseMatrix>>(400, 100, 80);
  TestRandomSystem<
      Group1130CudaDenseMatrix,
      Group1130CudaSparseMatrix,
      micm::CudaProcessSet<Group1130CudaDenseMatrix, Group1130CudaSparseMatrix>>(400, 100, 80);
}

TEST(CudaProcessSet, CudaMatrix)
{
  TestProcessSet<
      Group1CudaDenseMatrix,
      Group1CudaSparseMatrix,
      micm::CudaProcessSet<Group1CudaDenseMatrix, Group1CudaSparseMatrix>>();
  TestProcessSet<
      Group3CudaDenseMatrix,
      Group3CudaSparseMatrix,
      micm::CudaProcessSet<Group3CudaDenseMatrix, Group3CudaSparseMatrix>>();
  TestProcessSet<
      Group27CudaDenseMatrix,
      Group27CudaSparseMatrix,
      micm::CudaProcessSet<Group27CudaDenseMatrix, Group27CudaSparseMatrix>>();
  TestProcessSet<
      Group32CudaDenseMatrix,
      Group32CudaSparseMatrix,
      micm::CudaProcessSet<Group32CudaDenseMatrix, Group32CudaSparseMatrix>>();
  TestProcessSet<
      Group43CudaDenseMatrix,
      Group43CudaSparseMatrix,
      micm::CudaProcessSet<Group43CudaDenseMatrix, Group43CudaSparseMatrix>>();
  TestProcessSet<
      Group77CudaDenseMatrix,
      Group77CudaSparseMatrix,
      micm::CudaProcessSet<Group77CudaDenseMatrix, Group77CudaSparseMatrix>>();
  TestProcessSet<
      Group113CudaDenseMatrix,
      Group113CudaSparseMatrix,
      micm::CudaProcessSet<Group113CudaDenseMatrix, Group113CudaSparseMatrix>>();
  TestProcessSet<
      Group193CudaDenseMatrix,
      Group193CudaSparseMatrix,
      micm::CudaProcessSet<Group193CudaDenseMatrix, Group193CudaSparseMatrix>>();
  TestProcessSet<
      Group281CudaDenseMatrix,
      Group281CudaSparseMatrix,
      micm::CudaProcessSet<Group281CudaDenseMatrix, Group281CudaSparseMatrix>>();
  TestProcessSet<
      Group472CudaDenseMatrix,
      Group472CudaSparseMatrix,
      micm::CudaProcessSet<Group472CudaDenseMatrix, Group472CudaSparseMatrix>>();
  TestProcessSet<
      Group512CudaDenseMatrix,
      Group512CudaSparseMatrix,
      micm::CudaProcessSet<Group512CudaDenseMatrix, Group512CudaSparseMatrix>>();
  TestProcessSet<
      Group739CudaDenseMatrix,
      Group739CudaSparseMatrix,
      micm::CudaProcessSet<Group739CudaDenseMatrix, Group739CudaSparseMatrix>>();
  TestProcessSet<
      Group1130CudaDenseMatrix,
      Group1130CudaSparseMatrix,
      micm::CudaProcessSet<Group1130CudaDenseMatrix, Group1130CudaSparseMatrix>>();
}

TEST(CudaProcessSetAlgebraicVariables, CudaMatrix)
{
  testAlgebraicMasking<
      Group1CudaDenseMatrix,
      Group1CudaSparseMatrix,
      micm::CudaProcessSet<Group1CudaDenseMatrix, Group1CudaSparseMatrix>>();
  testAlgebraicMasking<
      Group3CudaDenseMatrix,
      Group3CudaSparseMatrix,
      micm::CudaProcessSet<Group3CudaDenseMatrix, Group3CudaSparseMatrix>>();
  testAlgebraicMasking<
      Group32CudaDenseMatrix,
      Group32CudaSparseMatrix,
      micm::CudaProcessSet<Group32CudaDenseMatrix, Group32CudaSparseMatrix>>();
  testAlgebraicMasking<
      Group43CudaDenseMatrix,
      Group43CudaSparseMatrix,
      micm::CudaProcessSet<Group43CudaDenseMatrix, Group43CudaSparseMatrix>>();
  testAlgebraicMasking<
      Group512CudaDenseMatrix,
      Group512CudaSparseMatrix,
      micm::CudaProcessSet<Group512CudaDenseMatrix, Group512CudaSparseMatrix>>();
}
