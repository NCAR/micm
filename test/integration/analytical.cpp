#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>
#include <utility>
#include <vector>

void writeCSV(const std::string& filename, const std::vector<std::vector<double>>& data)
{
  std::ofstream file(filename);
  if (file.is_open())
  {
    // Write column headers
    file << "time,A,B,C\n";

    // Write data rows
    for (size_t i = 0; i < data.size(); ++i)
    {
      file << i;  // Using i as the time value
      for (size_t j = 0; j < data[i].size(); ++j)
      {
        file << "," << data[i][j];
      }
      file << "\n";
    }
    file.close();
  }
  else
  {
    std::cerr << "Error opening file: " << filename << std::endl;
  }
}

using yields = std::pair<micm::Species, double>;

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

TEST(AnalyticalExamples, Troe)
{
  /*
   * A -k1-> B -k2-> C
   *
   * Copying the CAMP example: https://github.com/open-atmos/camp/blob/main/test/unit_rxn_data/test_rxn_troe.F90
   */
  constexpr size_t nsteps = 100;

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TroeRateConstant({ .k0_A_ = 4.0e-18 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::TroeRateConstant({ .k0_A_ = 1.2e-12,
                                                                 .k0_B_ = 167,
                                                                 .k0_C_ = 3,
                                                                 .kinf_A_ = 136,
                                                                 .kinf_B_ = 5,
                                                                 .kinf_C_ = 24,
                                                                 .Fc_ = 0.9,
                                                                 .N_ = 0.8 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1.0e6;

  // A->B reaction rate
  double k_0 = 4.0e-18;
  double k_inf = 1;
  double k1 = k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
        pow(0.6, 1.0 / (1.0 + (1.0 / 1.0) * pow(log10(k_0 * air_density / k_inf), 2)));
    
  // B->C reaction rate
  k_0 = 1.2e-12 * exp(3.0 / temperature) * pow(temperature / 300.0, 167.0);
  k_inf = 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  double k2 = k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
        pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * pow(log10(k_0 * air_density / k_inf), 2)));

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  std::vector<std::vector<double>> model_concentrations(nsteps, std::vector<double>(3));
  std::vector<std::vector<double>> analytical_concentrations(nsteps, std::vector<double>(3));

  model_concentrations[0] = { 1, 0, 0 };
  analytical_concentrations[0] = { 1, 0, 0 };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    // Model results
    auto result = solver.Solve(i_time, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    EXPECT_EQ(k1, state.rate_constants_.AsVector()[0]);
    EXPECT_EQ(k2, state.rate_constants_.AsVector()[1]);
    model_concentrations[i_time] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));
    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  writeCSV("analytical_concentrations.csv", analytical_concentrations);
  writeCSV("model_concentrations.csv", model_concentrations);

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    for (size_t j = 0; j < model_concentrations[i].size(); ++j)
    {
      EXPECT_DOUBLE_EQ(model_concentrations[i][j], analytical_concentrations[i][j])
          << "Arrays differ at index (" << i << ", " << j << ")";
    }
  }
}
