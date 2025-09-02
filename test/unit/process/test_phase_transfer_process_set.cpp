#include <micm/process/phase_transfer_process_set.hpp>

#include <gtest/gtest.h>

enum ForcingIndex
{
  IDX_NUMBER_CONCENTRATION = 1,
  IDX_EFFECTIVE_RADIUS,
};

class AerosolModel
{
 public:

  void AddForcingTerms(std::vector<std::vector<double>>& forcing)
  {
    // Set model specific parameters to update forcing terms
    forcing[0][ForcingIndex::IDX_NUMBER_CONCENTRATION] = 0.77777;
    forcing[0][ForcingIndex::IDX_EFFECTIVE_RADIUS] = 0.555555;

    std::cout << "[AerosolModel] Callback received from MICM" << std::endl;
  }
};

void TestPhaseTransferProcessSet() 
{
  micm::PhaseTransferProcessSet process_set;
  AerosolModel model;

  std::vector<std::vector<double>> forcing(3, std::vector<double>(3, 0.0));

  // Set the callback using lambda that calls the AerosolModel's function
  process_set.RegisterForcingCallback([&model](std::vector<std::vector<double>>& forcing){
      model.AddForcingTerms(forcing);
  });

  process_set.AddForcingTerms(forcing);

  EXPECT_EQ(forcing[0][0], 0.3333);
  EXPECT_EQ(forcing[0][ForcingIndex::IDX_NUMBER_CONCENTRATION], 0.77777);
  EXPECT_EQ(forcing[0][ForcingIndex::IDX_EFFECTIVE_RADIUS], 0.555555);
  EXPECT_EQ(forcing[1][0], 0.0);
  EXPECT_EQ(forcing[1][ForcingIndex::IDX_NUMBER_CONCENTRATION], 0.0);
  EXPECT_EQ(forcing[1][ForcingIndex::IDX_EFFECTIVE_RADIUS], 0.0);
}

TEST(PhaseTransferProcessSet, Forcing)
{
  TestPhaseTransferProcessSet();
}