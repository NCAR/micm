#include <micm/process/intraphase_process_builder.hpp>
#include <micm/system/species.hpp>
#include <micm/system/phase.hpp>

#include <gtest/gtest.h>

TEST(IntraPhaseProcessBuilder, DefaultConstructor){
  micm::IntraPhaseProcessBuilder builder{};
}

TEST(IntraPhaseProcessBuilder, CanChainBuilder){
  micm::IntraPhaseProcessBuilder builder{};

  builder
    .For(micm::Phase())
    .With(micm::Species());
}