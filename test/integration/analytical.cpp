#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>
#include <utility>
#include <vector>

constexpr size_t nsteps = 1000;

void writeCSV(
    const std::string& filename,
    const std::vector<std::string>& header,
    const std::vector<std::vector<double>>& data)
{
  std::ofstream file(filename);
  if (file.is_open())
  {
    // Write column headers
    for (size_t i = 0; i < header.size(); ++i)
    {
      file << header[i];
      if (i < header.size() - 1)
      {
        file << ",";
      }
    }
    file << "\n";

    // Write data rows
    for (size_t i = 0; i < data.size(); ++i)
    {
      file << i << ",";
      for (size_t j = 0; j < data[i].size(); ++j)
      {
        file << data[i][j];
        if (j < data[i].size() - 1)
        {
          file << ",";
        }
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
   * A -> B, k1
   * B -> C, k2
   *
   * Copying the CAMP example: https://github.com/open-atmos/camp/blob/main/test/unit_rxn_data/test_rxn_troe.F90
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TroeRateConstant({ .k0_A_ = 4.0e-10 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::TroeRateConstant({ .k0_A_ = 1.2e-3,
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
  double air_density = 1e6;

  // A->B reaction rate
  double k_0 = 4.0e-10;
  double k_inf = 1;
  double k1 = k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
              pow(0.6, 1.0 / (1.0 + (1.0 / 1.0) * pow(log10(k_0 * air_density / k_inf), 2)));

  // B->C reaction rate
  k_0 = 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 167.0);
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
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    EXPECT_NEAR(k1, state.rate_constants_.AsVector()[0], 1e-8);
    EXPECT_NEAR(k2, state.rate_constants_.AsVector()[1], 1e-8);
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations);
  writeCSV("model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a], analytical_concentrations[i][0], 1e-8)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-8)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-8)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  /*
   * A1 -> B, k1
   * A2 -> B, k1
   * A1 -> A2, k3 >>> k1
   * A2 -> A1, k4 >>> k1
   * B -> C, k2
   *
   */

  auto a1 = micm::Species("A1");
  auto a2 = micm::Species("A2");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a1, a2, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TroeRateConstant({ .k0_A_ = 4.0e-10 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TroeRateConstant({ .k0_A_ = 4.0e-10 }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::TroeRateConstant({ .k0_A_ = 1.2e-3,
                                                                 .k0_B_ = 167,
                                                                 .k0_C_ = 3,
                                                                 .kinf_A_ = 136,
                                                                 .kinf_B_ = 5,
                                                                 .kinf_C_ = 24,
                                                                 .Fc_ = 0.9,
                                                                 .N_ = 0.8 }))
                         .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(a2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .phase(gas_phase);

  micm::Process r5 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(a1, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2, r3, r4, r5 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k_0 = 4.0e-10;
  double k_inf = 1;
  double k1 = k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
              pow(0.6, 1.0 / (1.0 + (1.0 / 1.0) * pow(log10(k_0 * air_density / k_inf), 2)));

  // B->C reaction rate
  k_0 = 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 167.0);
  k_inf = 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  double k2 = k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
              pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * pow(log10(k_0 * air_density / k_inf), 2)));

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  std::vector<std::vector<double>> model_concentrations(nsteps, std::vector<double>(4));
  std::vector<std::vector<double>> analytical_concentrations(nsteps, std::vector<double>(3));

  state.SetConcentration(a1, 0.5);
  state.SetConcentration(a2, 0.5);
  state.SetConcentration(b, 0.0);
  state.SetConcentration(c, 0.0);

  model_concentrations[0] = state.variables_[0];
  analytical_concentrations[0] = { 1, 0, 0 };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    // Model results
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  auto header = state.variable_names_;
  header.insert(header.begin(), "time");
  writeCSV("stiff_model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a1 = map.at("A1");
  size_t _a2 = map.at("A2");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a1] + model_concentrations[i][_a2], analytical_concentrations[i][0], 1e-4);
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-4);
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-4);
  }
}

TEST(AnalyticalExamples, Photolysis)
{
  /*
   * A -> B, k1
   * B -> C, k2
   *
   * Copying the CAMP example: https://github.com/open-atmos/camp/blob/main/test/unit_rxn_data/test_rxn_photolysis.F90
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "photoA" }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "photoB" }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k1 = 2e-3;

  // B->C reaction rate
  double k2 = 3e-3;

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  state.SetCustomRateParameter("photoA", k1);
  state.SetCustomRateParameter("photoB", k2);

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
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    EXPECT_NEAR(k1, state.rate_constants_.AsVector()[0], 1e-8);
    EXPECT_NEAR(k2, state.rate_constants_.AsVector()[1], 1e-8);
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations);
  writeCSV("model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a], analytical_concentrations[i][0], 1e-8)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-8)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-8)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  /*
   * A1 -> B, k1
   * A2 -> B, k1
   * A1 -> A2, k3 >>> k1
   * A2 -> A1, k4 >>> k1
   * B -> C, k2
   *
   */

  auto a1 = micm::Species("A1");
  auto a2 = micm::Species("A2");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a1, a2, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "photoA1B" }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "photoA2B" }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "photoB" }))
                         .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(a2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .phase(gas_phase);

  micm::Process r5 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(a1, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2, r3, r4, r5 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k1 = 2e-3;

  // B->C reaction rate
  double k2 = 3e-3;

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  state.SetCustomRateParameter("photoA1B", k1);
  state.SetCustomRateParameter("photoA2B", k1);
  state.SetCustomRateParameter("photoB", k2);

  std::vector<std::vector<double>> model_concentrations(nsteps, std::vector<double>(3));
  std::vector<std::vector<double>> analytical_concentrations(nsteps, std::vector<double>(3));

  state.SetConcentration(a1, 0.5);
  state.SetConcentration(a2, 0.5);
  state.SetConcentration(b, 0.0);
  state.SetConcentration(c, 0.0);

  model_concentrations[0] = state.variables_[0];

  analytical_concentrations[0] = { 1, 0, 0 };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    // Model results
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations);
  writeCSV("model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a1 = map.at("A1");
  size_t _a2 = map.at("A2");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a1] + model_concentrations[i][_a2], analytical_concentrations[i][0], 1e-4);
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-4);
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-4);
  }
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  /*
   * A -> B, k1
   * B -> C, k2
   *
   * Copying the CAMP example:
   * https://github.com/open-atmos/camp/blob/main/test/unit_rxn_data/test_rxn_ternary_chemical_activation.F90
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 4.0e-10, .kinf_A_ = 1 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2e-3,
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
  double air_density = 1e6;

  // A->B reaction rate
  double k_0 = 4.0e-10;
  double k_inf = 1;
  double k1 = k_0 / (1.0 + k_0 * air_density / k_inf) *
              pow(0.6, 1.0 / (1.0 + (1.0 / 1.0) * pow(log10(k_0 * air_density / k_inf), 2)));

  // B->C reaction rate
  k_0 = 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 167.0);
  k_inf = 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  double k2 = k_0 / (1.0 + k_0 * air_density / k_inf) *
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
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    EXPECT_NEAR(k1, state.rate_constants_.AsVector()[0], 1e-8);
    EXPECT_NEAR(k2, state.rate_constants_.AsVector()[1], 1e-8);
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations);
  writeCSV("model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a], analytical_concentrations[i][0], 1e-8)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-8)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-8)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  /*
   * A1 -> B, k1
   * A2 -> B, k1
   * A1 -> A2, k3 >>> k1
   * A2 -> A1, k4 >>> k1
   * B -> C, k2
   *
   */

  auto a1 = micm::Species("A1");
  auto a2 = micm::Species("A2");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a1, a2, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 4.0e-10, .kinf_A_ = 1 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 4.0e-10, .kinf_A_ = 1 }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2e-3,
                                                                                      .k0_B_ = 167,
                                                                                      .k0_C_ = 3,
                                                                                      .kinf_A_ = 136,
                                                                                      .kinf_B_ = 5,
                                                                                      .kinf_C_ = 24,
                                                                                      .Fc_ = 0.9,
                                                                                      .N_ = 0.8 }))
                         .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(a2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .phase(gas_phase);

  micm::Process r5 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(a1, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2, r3, r4, r5 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k_0 = 4.0e-10;
  double k_inf = 1;
  double k1 = k_0 / (1.0 + k_0 * air_density / k_inf) *
              pow(0.6, 1.0 / (1.0 + (1.0 / 1.0) * pow(log10(k_0 * air_density / k_inf), 2)));

  // B->C reaction rate
  k_0 = 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 167.0);
  k_inf = 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  double k2 = k_0 / (1.0 + k_0 * air_density / k_inf) *
              pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * pow(log10(k_0 * air_density / k_inf), 2)));

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  std::vector<std::vector<double>> model_concentrations(nsteps, std::vector<double>(4));
  std::vector<std::vector<double>> analytical_concentrations(nsteps, std::vector<double>(3));

  state.SetConcentration(a1, 0.5);
  state.SetConcentration(a2, 0.5);
  state.SetConcentration(b, 0.0);
  state.SetConcentration(c, 0.0);

  model_concentrations[0] = state.variables_[0];
  analytical_concentrations[0] = { 1, 0, 0 };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    // Model results
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  auto header = state.variable_names_;
  header.insert(header.begin(), "time");
  writeCSV("stiff_model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a1 = map.at("A1");
  size_t _a2 = map.at("A2");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a1] + model_concentrations[i][_a2], analytical_concentrations[i][0], 1e-4);
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-4);
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-4);
  }
}

TEST(AnalyticalExamples, Tunneling)
{
  /*
   * A -> B, k1
   * B -> C, k2
   *
   * Copying the CAMP example:
   * https://github.com/open-atmos/camp/blob/main/test/unit_rxn_data/test_rxn_wennberg_tunneling.F90
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TunnelingRateConstant({ .A_ = 4.0e-3 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::TunnelingRateConstant({ .A_ = 1.2e-4, .B_ = 167, .C_ = 1.0e8 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k1 = 4.0e-3;

  // B->C reaction rate
  double k2 = 1.2e-4 * exp(-167 / temperature + 1.0e8 / pow(temperature, 3));

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
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    EXPECT_NEAR(k1, state.rate_constants_.AsVector()[0], 1e-8);
    EXPECT_NEAR(k2, state.rate_constants_.AsVector()[1], 1e-8);
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations);
  writeCSV("model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a], analytical_concentrations[i][0], 1e-8)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-8)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-8)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  /*
   * A1 -> B, k1
   * A2 -> B, k1
   * A1 -> A2, k3 >>> k1
   * A2 -> A1, k4 >>> k1
   * B -> C, k2
   *
   */

  auto a1 = micm::Species("A1");
  auto a2 = micm::Species("A2");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a1, a2, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TunnelingRateConstant({ .A_ = 4.0e-3 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::TunnelingRateConstant({ .A_ = 4.0e-3 }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::TunnelingRateConstant({ .A_ = 1.2e-4, .B_ = 167, .C_ = 1.0e8 }))
                         .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(a2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .phase(gas_phase);

  micm::Process r5 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(a1, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2, r3, r4, r5 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k1 = 4.0e-3;

  // B->C reaction rate
  double k2 = 1.2e-4 * exp(-167 / temperature + 1.0e8 / pow(temperature, 3));

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  std::vector<std::vector<double>> model_concentrations(nsteps, std::vector<double>(4));
  std::vector<std::vector<double>> analytical_concentrations(nsteps, std::vector<double>(3));

  state.SetConcentration(a1, 0.5);
  state.SetConcentration(a2, 0.5);
  state.SetConcentration(b, 0.0);
  state.SetConcentration(c, 0.0);

  model_concentrations[0] = state.variables_[0];
  analytical_concentrations[0] = { 1, 0, 0 };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    // Model results
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  auto header = state.variable_names_;
  header.insert(header.begin(), "time");
  writeCSV("stiff_model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a1 = map.at("A1");
  size_t _a2 = map.at("A2");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a1] + model_concentrations[i][_a2], analytical_concentrations[i][0], 1e-4)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-4)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-4)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, Arrhenius)
{
  /*
   * A -> B, k1
   * B -> C, k2
   *
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e-3, .C_ = 50 }))
                         .phase(gas_phase);

  micm::Process r2 =
      micm::Process::create()
          .reactants({ b })
          .products({ yields(c, 1) })
          .rate_constant(micm::ArrheniusRateConstant({ .A_ = 1.2e-4, .B_ = 167, .C_ = 75, .D_ = 50, .E_ = 0.5 }))
          .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k1 = 4.0e-3 * std::exp(50 / temperature);

  // B->C reaction rate
  double k2 = 1.2e-4 * std::exp(75 / temperature) * pow(temperature / 50, 167) * (1.0 + 0.5 * pressure);

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
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    EXPECT_NEAR(k1, state.rate_constants_.AsVector()[0], 1e-8);
    EXPECT_NEAR(k2, state.rate_constants_.AsVector()[1], 1e-8);
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations);
  writeCSV("model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a], analytical_concentrations[i][0], 1e-8)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-8)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-8)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  /*
   * A1 -> B, k1
   * A2 -> B, k1
   * A1 -> A2, k3 >>> k1
   * A2 -> A1, k4 >>> k1
   * B -> C, k2
   *
   */

  auto a1 = micm::Species("A1");
  auto a2 = micm::Species("A2");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a1, a2, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e-3, .C_ = 50 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e-3, .C_ = 50 }))
                         .phase(gas_phase);

  micm::Process r3 =
      micm::Process::create()
          .reactants({ b })
          .products({ yields(c, 1) })
          .rate_constant(micm::ArrheniusRateConstant({ .A_ = 1.2e-4, .B_ = 167, .C_ = 75, .D_ = 50, .E_ = 0.5 }))
          .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(a2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .phase(gas_phase);

  micm::Process r5 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(a1, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2, r3, r4, r5 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double k1 = 4.0e-3 * std::exp(50 / temperature);

  // B->C reaction rate
  double k2 = 1.2e-4 * std::exp(75 / temperature) * pow(temperature / 50, 167) * (1.0 + 0.5 * pressure);

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  std::vector<std::vector<double>> model_concentrations(nsteps, std::vector<double>(4));
  std::vector<std::vector<double>> analytical_concentrations(nsteps, std::vector<double>(3));

  state.SetConcentration(a1, 0.5);
  state.SetConcentration(a2, 0.5);
  state.SetConcentration(b, 0.0);
  state.SetConcentration(c, 0.0);

  model_concentrations[0] = state.variables_[0];
  analytical_concentrations[0] = { 1, 0, 0 };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    // Model results
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  auto header = state.variable_names_;
  header.insert(header.begin(), "time");
  writeCSV("stiff_model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a1 = map.at("A1");
  size_t _a2 = map.at("A2");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a1] + model_concentrations[i][_a2], analytical_concentrations[i][0], 1e-4)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-4)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-4)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, Branched)
{
  /*
   * A -> B, k1
   * B -> C, k2
   *
   * Copying the CAMP example: https://github.com/open-atmos/camp/blob/main/test/unit_rxn_data/test_rxn_wennberg_no_ro2.F90
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 =
      micm::Process::create()
          .reactants({ a })
          .products({ yields(b, 1) })
          .rate_constant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy,
                                                      .X_ = 1.2,
                                                      .Y_ = 204.3,
                                                      .a0_ = 1.0e-3,
                                                      .n_ = 2 }))
          .phase(gas_phase);

  micm::Process r2 =
      micm::Process::create()
          .reactants({ b })
          .products({ yields(c, 1) })
          .rate_constant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Nitrate,
                                                      .X_ = 1.2,
                                                      .Y_ = 204.3,
                                                      .a0_ = 1.0e-3,
                                                      .n_ = 2 }))
          .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double air_dens_n_cm3 = air_density * AVOGADRO_CONSTANT * 1.0e-6;
  double a_ = 2.0e-22 * std::exp(2) * 2.45e19;
  double b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
  double z = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
  double A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

  double k1 = 1.2 * std::exp(-204.3 / temperature) * (z / (z + A));

  // B->C reaction rate
  a_ = 2.0e-22 * std::exp(2) * 2.45e19;
  b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
  z = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
  A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

  double k2 = 1.2 * std::exp(-204.3 / temperature) * (A / (z + A));

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
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    EXPECT_NEAR(k1, state.rate_constants_.AsVector()[0], 1e-8);
    EXPECT_NEAR(k2, state.rate_constants_.AsVector()[1], 1e-8);
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations);
  writeCSV("model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a], analytical_concentrations[i][0], 1e-3)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-3)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-3)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  /*
   * A1 -> B, k1
   * A2 -> B, k1
   * A1 -> A2, k3 >>> k1
   * A2 -> A1, k4 >>> k1
   * B -> C, k2
   *
   */

  auto a1 = micm::Species("A1");
  auto a2 = micm::Species("A2");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a1, a2, b, c } };

  micm::Process r1 =
      micm::Process::create()
          .reactants({ a1 })
          .products({ yields(b, 1) })
          .rate_constant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy,
                                                      .X_ = 1.2,
                                                      .Y_ = 204.3,
                                                      .a0_ = 1.0e-3,
                                                      .n_ = 2 }))
          .phase(gas_phase);

  micm::Process r2 =
      micm::Process::create()
          .reactants({ a2 })
          .products({ yields(b, 1) })
          .rate_constant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy,
                                                      .X_ = 1.2,
                                                      .Y_ = 204.3,
                                                      .a0_ = 1.0e-3,
                                                      .n_ = 2 }))
          .phase(gas_phase);

  micm::Process r3 =
      micm::Process::create()
          .reactants({ b })
          .products({ yields(c, 1) })
          .rate_constant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Nitrate,
                                                      .X_ = 1.2,
                                                      .Y_ = 204.3,
                                                      .a0_ = 1.0e-3,
                                                      .n_ = 2 }))
          .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ a1 })
                         .products({ yields(a2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .phase(gas_phase);

  micm::Process r5 = micm::Process::create()
                         .reactants({ a2 })
                         .products({ yields(a1, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2, r3, r4, r5 },
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  // A->B reaction rate
  double air_dens_n_cm3 = air_density * AVOGADRO_CONSTANT * 1.0e-6;
  double a_ = 2.0e-22 * std::exp(2) * 2.45e19;
  double b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
  double z = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
  double A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

  double k1 = 1.2 * std::exp(-204.3 / temperature) * (z / (z + A));

  // B->C reaction rate
  a_ = 2.0e-22 * std::exp(2) * 2.45e19;
  b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
  z = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
  A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

  double k2 = 1.2 * std::exp(-204.3 / temperature) * (A / (z + A));

  double time_step = 1.0;
  micm::State<micm::Matrix> state = solver.GetState();

  std::vector<std::vector<double>> model_concentrations(nsteps, std::vector<double>(4));
  std::vector<std::vector<double>> analytical_concentrations(nsteps, std::vector<double>(3));

  state.SetConcentration(a1, 0.5);
  state.SetConcentration(a2, 0.5);
  state.SetConcentration(b, 0.0);
  state.SetConcentration(c, 0.0);

  model_concentrations[0] = state.variables_[0];
  analytical_concentrations[0] = { 1, 0, 0 };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    // Model results
    auto result = solver.Solve(time_step, state);
    EXPECT_EQ(result.state_, (micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>::SolverState::Converged));
    model_concentrations[i_time] = result.result_.AsVector();
    state.variables_[0] = result.result_.AsVector();

    // Analytical results
    double time = i_time * time_step;

    double initial_A = analytical_concentrations[0][idx_A];
    analytical_concentrations[i_time][idx_A] = initial_A * std::exp(-(k1)*time);
    analytical_concentrations[i_time][idx_B] = initial_A * (k1 / (k2 - k1)) * (std::exp(-k1 * time) - std::exp(-k2 * time));

    analytical_concentrations[i_time][idx_C] =
        initial_A * (1.0 + (k1 * std::exp(-k2 * time) - k2 * std::exp(-k1 * time)) / (k2 - k1));
  }

  auto header = state.variable_names_;
  header.insert(header.begin(), "time");
  writeCSV("stiff_model_concentrations.csv", header, model_concentrations);

  auto map = state.variable_map_;

  size_t _a1 = map.at("A1");
  size_t _a2 = map.at("A2");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_a1] + model_concentrations[i][_a2], analytical_concentrations[i][0], 1e-4)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_b], analytical_concentrations[i][1], 1e-4)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_c], analytical_concentrations[i][2], 1e-4)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}