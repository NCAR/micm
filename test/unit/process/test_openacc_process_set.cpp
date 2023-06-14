#include <gtest/gtest.h>

#include <micm/process/process.hpp>
#include <micm/process/process_set_openacc.hpp>

TEST(ProcessSet, Deriv)
{
  micm::openacc::deriv();
}