#include <micm/process/intraphase_process_builder.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/external_rate_constant.hpp>
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
    .With(micm::Species())
    .Producing(micm::Species())
    .WithYield(2l)
    .WithRateConstant(micm::ArrheniusRateConstant())
    .WithRateConstant(micm::TroeRateConstant())
    .WithRateConstant(micm::ExternalRateConstant())
    .Build();
}