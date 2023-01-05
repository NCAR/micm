#include <micm/process/intraphase_process_builder.hpp>
#include <micm/system/species.hpp>
#include <micm/system/phase.hpp>

#include <gtest/gtest.h>

TEST(IntraPhaseProcessBuilder, DefaultConstructor){
  micm::IntraPhaseProcessBuilder<double> builder{};
}

TEST(IntraPhaseProcessBuilder, CanChainBuilder){
  micm::IntraPhaseProcessBuilder<double> builder{};

  builder
    .For(micm::Phase<double>())
    .With(micm::Species<double>());
}