#include <micm/process/intraphase_process.hpp>

#include <gtest/gtest.h>

TEST(IntraPhaseProcess, DefaultConstructor){
  micm::IntraPhaseProcess<double> process{};
}