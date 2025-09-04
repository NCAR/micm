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
  testRandomSystem<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group3CudaDenseMatrix, Group3CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group27CudaDenseMatrix, Group27CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group32CudaDenseMatrix, Group32CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group43CudaDenseMatrix, Group43CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group77CudaDenseMatrix, Group77CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group113CudaDenseMatrix, Group113CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group193CudaDenseMatrix, Group193CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group281CudaDenseMatrix, Group281CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group472CudaDenseMatrix, Group472CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group512CudaDenseMatrix, Group512CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group739CudaDenseMatrix, Group739CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
  testRandomSystem<Group1130CudaDenseMatrix, Group1130CudaSparseMatrix, micm::CudaProcessSet>(400, 100, 80);
}

TEST(CudaProcessSet, CudaMatrix)
{
  testProcessSet<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group3CudaDenseMatrix, Group3CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group27CudaDenseMatrix, Group27CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group32CudaDenseMatrix, Group32CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group43CudaDenseMatrix, Group43CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group77CudaDenseMatrix, Group77CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group193CudaDenseMatrix, Group193CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group281CudaDenseMatrix, Group281CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group472CudaDenseMatrix, Group472CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group512CudaDenseMatrix, Group512CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group113CudaDenseMatrix, Group113CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group739CudaDenseMatrix, Group739CudaSparseMatrix, micm::CudaProcessSet>();
  testProcessSet<Group1130CudaDenseMatrix, Group1130CudaSparseMatrix, micm::CudaProcessSet>();
}
