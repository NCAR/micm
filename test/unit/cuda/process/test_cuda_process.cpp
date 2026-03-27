// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "../../process/test_process_policy.hpp"

#include <micm/cuda/util/cuda_dense_matrix.hpp>

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

TEST(CudaProcess, ArrheniusRateConstants)
{
  testCudaProcessRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, MixedRateConstants)
{
  testCudaProcessMixedRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessMixedRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessMixedRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessMixedRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessMixedRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessMixedRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessMixedRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessMixedRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessMixedRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessMixedRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessMixedRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessMixedRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessMixedRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, TroeRateConstants)
{
  testCudaProcessTroeRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessTroeRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessTroeRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessTroeRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessTroeRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessTroeRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessTroeRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessTroeRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessTroeRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessTroeRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessTroeRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessTroeRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessTroeRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, TunnelingRateConstants)
{
  testCudaProcessTunnelingRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessTunnelingRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessTunnelingRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessTunnelingRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessTunnelingRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessTunnelingRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessTunnelingRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessTunnelingRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessTunnelingRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessTunnelingRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessTunnelingRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessTunnelingRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessTunnelingRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, BranchedAlkoxyRateConstants)
{
  testCudaProcessBranchedAlkoxyRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessBranchedAlkoxyRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessBranchedAlkoxyRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessBranchedAlkoxyRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessBranchedAlkoxyRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessBranchedAlkoxyRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessBranchedAlkoxyRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessBranchedAlkoxyRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessBranchedAlkoxyRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessBranchedAlkoxyRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessBranchedAlkoxyRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessBranchedAlkoxyRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessBranchedAlkoxyRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, BranchedNitrateRateConstants)
{
  testCudaProcessBranchedNitrateRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessBranchedNitrateRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessBranchedNitrateRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessBranchedNitrateRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessBranchedNitrateRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessBranchedNitrateRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessBranchedNitrateRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessBranchedNitrateRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessBranchedNitrateRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessBranchedNitrateRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessBranchedNitrateRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessBranchedNitrateRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessBranchedNitrateRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, TernaryChemicalActivationRateConstants)
{
  testCudaProcessTernaryChemicalActivationRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessTernaryChemicalActivationRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessTernaryChemicalActivationRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessTernaryChemicalActivationRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessTernaryChemicalActivationRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessTernaryChemicalActivationRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessTernaryChemicalActivationRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessTernaryChemicalActivationRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessTernaryChemicalActivationRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessTernaryChemicalActivationRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessTernaryChemicalActivationRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessTernaryChemicalActivationRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessTernaryChemicalActivationRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, ReversibleRateConstants)
{
  testCudaProcessReversibleRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessReversibleRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessReversibleRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessReversibleRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessReversibleRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessReversibleRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessReversibleRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessReversibleRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessReversibleRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessReversibleRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessReversibleRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessReversibleRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessReversibleRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, SurfaceRateConstants)
{
  testCudaProcessSurfaceRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessSurfaceRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessSurfaceRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessSurfaceRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessSurfaceRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessSurfaceRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessSurfaceRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessSurfaceRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessSurfaceRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessSurfaceRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessSurfaceRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessSurfaceRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessSurfaceRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, UserDefinedRateConstants)
{
  testCudaProcessUserDefinedRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessUserDefinedRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessUserDefinedRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessUserDefinedRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessUserDefinedRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessUserDefinedRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessUserDefinedRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessUserDefinedRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessUserDefinedRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessUserDefinedRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessUserDefinedRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessUserDefinedRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessUserDefinedRateConstants<Group1130CudaDenseMatrix>(1200);
}

TEST(CudaProcess, AllRateConstants)
{
  testCudaProcessAllRateConstants<Group1CudaDenseMatrix>(5);
  testCudaProcessAllRateConstants<Group3CudaDenseMatrix>(5);
  testCudaProcessAllRateConstants<Group27CudaDenseMatrix>(50);
  testCudaProcessAllRateConstants<Group32CudaDenseMatrix>(50);
  testCudaProcessAllRateConstants<Group43CudaDenseMatrix>(50);
  testCudaProcessAllRateConstants<Group77CudaDenseMatrix>(100);
  testCudaProcessAllRateConstants<Group113CudaDenseMatrix>(200);
  testCudaProcessAllRateConstants<Group193CudaDenseMatrix>(200);
  testCudaProcessAllRateConstants<Group281CudaDenseMatrix>(400);
  testCudaProcessAllRateConstants<Group472CudaDenseMatrix>(500);
  testCudaProcessAllRateConstants<Group512CudaDenseMatrix>(600);
  testCudaProcessAllRateConstants<Group739CudaDenseMatrix>(800);
  testCudaProcessAllRateConstants<Group1130CudaDenseMatrix>(1200);
}
