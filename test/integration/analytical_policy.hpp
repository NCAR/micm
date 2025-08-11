#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

constexpr size_t nsteps = 1000;
constexpr size_t NUM_CELLS = 3;

///////////////////////////
// Common test functions //
///////////////////////////

double relative_error(double a, double b)
{
  return abs(a - b) / abs(b);
}

double relative_difference(double a, double b)
{
  return abs(a - b) / ((a + b) / 2);
}

double combined_error(double a, double b, double abs_tol)
{
  return abs(a - b) * 2 / (abs(a) + abs(b) + abs_tol);
}

void writeCSV(
    const std::string& filename,
    const std::vector<std::string>& header,
    const std::vector<std::vector<double>>& data,
    const std::vector<double>& times)
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
      file << times[i] << ",";
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

void writeCSV2D(
    const std::string& filename,
    const std::vector<std::string>& header,
    const std::vector<std::vector<std::vector<double>>>& data,
    const std::vector<double>& times)
{
  std::ofstream file(filename);
  if (file.is_open())
  {
    // Write column headers
    for (std::size_t i = 0; i < header.size(); ++i)
    {
      file << header[i];
      if (i < header.size() - 1)
      {
        file << ",";
      }
    }
    file << "\n";

    // Write data rows
    for (std::size_t i_time = 0; i_time < times.size(); ++i_time)
    {
      for (std::size_t i_cell = 0; i_cell < data[i_time].size(); ++i_cell)
      {
        file << times[i_time] << "," << i_cell << ",";
        for (size_t j = 0; j < data[i_time][i_cell].size(); ++j)
        {
          file << data[i_time][i_cell][j];
          if (j < data[i_time][i_cell].size() - 1)
          {
            file << ",";
          }
        }
        file << "\n";
      }
    }
    file.close();
  }
  else
  {
    std::cerr << "Error opening file: " << filename << std::endl;
  }
}

double calculate_air_density_mol_m3(double pressure, double temperature)
{
  return pressure / (micm::constants::GAS_CONSTANT * temperature);
}

using SparseMatrixTest = micm::SparseMatrix<double>;

// Test the analytical solution for a simple A -k1-> B -k2-> C system
template<class BuilderPolicy>
void test_simple_system(
    const std::string& test_label,
    BuilderPolicy builder,
    double absolute_tolerances,
    std::function<double(double temperature, double pressure, double air_density)> calculate_k1,
    std::function<double(double temperature, double pressure, double air_density)> calculate_k2,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve,
    std::unordered_map<std::string, std::vector<double>> custom_parameters = {})
{
  auto solver = builder.Build();

  std::vector<double> temperatures = { 272.5, 254.7, 312.6 };
  std::vector<double> pressures = { 101253.3, 100672.5, 101319.8 };

  std::vector<double> k1(NUM_CELLS), k2(NUM_CELLS);

  for (int i = 0; i < NUM_CELLS; ++i)
  {
    double temperature = temperatures[i];
    double pressure = pressures[i];
    double air_density = calculate_air_density_mol_m3(pressure, temperature);
    k1[i] = calculate_k1(temperature, pressure, air_density);
    k2[i] = calculate_k2(temperature, pressure, air_density);
  }

  double time_step = 1.0;
  auto state = solver.GetState(NUM_CELLS);
  auto map = state.variable_map_;

  state.SetCustomRateParameters(custom_parameters);

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  std::vector<std::vector<std::vector<double>>> model_concentrations(
      nsteps, std::vector<std::vector<double>>(NUM_CELLS, std::vector<double>(3)));
  std::vector<std::vector<std::vector<double>>> analytical_concentrations(
      nsteps, std::vector<std::vector<double>>(NUM_CELLS, std::vector<double>(3)));

  for (int i = 0; i < NUM_CELLS; ++i)
  {
    model_concentrations[0][i][idx_A] = 1.0 - (double)i / (double)NUM_CELLS;
    model_concentrations[0][i][idx_B] = 0.0;
    model_concentrations[0][i][idx_C] = 0.0;
    analytical_concentrations[0][i] = model_concentrations[0][i];

    state.variables_[i][_a] = model_concentrations[0][i][idx_A];
    state.variables_[i][_b] = model_concentrations[0][i][idx_B];
    state.variables_[i][_c] = model_concentrations[0][i][idx_C];
    state.conditions_[i].temperature_ = temperatures[i];
    state.conditions_[i].pressure_ = pressures[i];
    state.conditions_[i].air_density_ = calculate_air_density_mol_m3(pressures[i], temperatures[i]);
  }

  std::vector<double> times;
  times.push_back(0);
  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    solver.CalculateRateConstants(state);
    prepare_for_solve(state);
    // Model results
    auto result = solver.Solve(time_step, state);
    postpare_for_solve(state);
    EXPECT_EQ(result.state_, (micm::SolverState::Converged));
    for (std::size_t i = 0; i < NUM_CELLS; ++i)
    {
      EXPECT_NEAR(k1[i], state.rate_constants_[i][0], absolute_tolerances);
      EXPECT_NEAR(k2[i], state.rate_constants_[i][1], absolute_tolerances);
      model_concentrations[i_time][i][idx_A] = state.variables_[i][_a];
      model_concentrations[i_time][i][idx_B] = state.variables_[i][_b];
      model_concentrations[i_time][i][idx_C] = state.variables_[i][_c];
    }

    // Analytical results
    double time = i_time * time_step;
    times.push_back(time);

    for (std::size_t i = 0; i < NUM_CELLS; ++i)
    {
      double initial_A = analytical_concentrations[0][i][idx_A];
      analytical_concentrations[i_time][i][idx_A] = initial_A * std::exp(-(k1[i]) * time);
      analytical_concentrations[i_time][i][idx_B] =
          initial_A * (k1[i] / (k2[i] - k1[i])) * (std::exp(-k1[i] * time) - std::exp(-k2[i] * time));

      analytical_concentrations[i_time][i][idx_C] =
          initial_A * (1.0 + (k1[i] * std::exp(-k2[i] * time) - k2[i] * std::exp(-k1[i] * time)) / (k2[i] - k1[i]));
    }
  }

  std::vector<std::string> header = { "time", "cell", "A", "B", "C" };
  writeCSV2D(test_label + "_analytical_concentrations.csv", header, analytical_concentrations, times);
  writeCSV2D(test_label + "_model_concentrations.csv", header, model_concentrations, times);

  for (std::size_t i_time = 1; i_time < model_concentrations.size(); ++i_time)
  {
    for (std::size_t i_cell = 0; i_cell < model_concentrations[i_time].size(); ++i_cell)
    {
      EXPECT_NEAR(
          model_concentrations[i_time][i_cell][idx_A], analytical_concentrations[i_time][i_cell][idx_A], absolute_tolerances)
          << "Arrays differ at index (" << i_time << ", " << i_cell << ", " << 0 << ") for " << test_label;
      EXPECT_NEAR(
          model_concentrations[i_time][i_cell][idx_B], analytical_concentrations[i_time][i_cell][idx_B], absolute_tolerances)
          << "Arrays differ at index (" << i_time << ", " << i_cell << ", " << 1 << ") for " << test_label;
      EXPECT_NEAR(
          model_concentrations[i_time][i_cell][idx_C], analytical_concentrations[i_time][i_cell][idx_C], absolute_tolerances)
          << "Arrays differ at index (" << i_time << ", " << i_cell << ", " << 2 << ") for " << test_label;
    }
  }
}

// Test the analytical solution for a simple stiff A1<-fast->A2 -k1-> B -k2-> C system
template<class BuilderPolicy>
void test_simple_stiff_system(
    const std::string& test_label,
    BuilderPolicy builder,
    double absolute_tolerances,
    std::function<double(double temperature, double pressure, double air_density)> calculate_k1,
    std::function<double(double temperature, double pressure, double air_density)> calculate_k2,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve,
    std::unordered_map<std::string, std::vector<double>> custom_parameters = {})
{
  auto solver = builder.Build();

  std::vector<double> temperatures = { 272.5, 254.7, 312.6 };
  std::vector<double> pressures = { 101253.3, 100672.5, 101319.8 };

  std::vector<double> k1(NUM_CELLS), k2(NUM_CELLS);

  for (int i = 0; i < NUM_CELLS; ++i)
  {
    double temperature = temperatures[i];
    double pressure = pressures[i];
    double air_density = calculate_air_density_mol_m3(pressure, temperature);
    k1[i] = calculate_k1(temperature, pressure, air_density);
    k2[i] = calculate_k2(temperature, pressure, air_density);
  }

  double time_step = 1.0;
  auto state = solver.GetState(NUM_CELLS);
  auto map = state.variable_map_;

  state.SetCustomRateParameters(custom_parameters);

  size_t idx_A = 0, idx_B = 1, idx_C = 2;

  size_t _a1 = map.at("A1");
  size_t _a2 = map.at("A2");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  std::vector<std::vector<std::vector<double>>> model_concentrations(
      nsteps, std::vector<std::vector<double>>(NUM_CELLS, std::vector<double>(3)));
  std::vector<std::vector<std::vector<double>>> analytical_concentrations(
      nsteps, std::vector<std::vector<double>>(NUM_CELLS, std::vector<double>(3)));

  for (int i = 0; i < NUM_CELLS; ++i)
  {
    model_concentrations[0][i][idx_A] = 1.0 - (double)i / (double)NUM_CELLS;
    model_concentrations[0][i][idx_B] = 0.0;
    model_concentrations[0][i][idx_C] = 0.0;
    analytical_concentrations[0][i] = model_concentrations[0][i];

    state.variables_[i][_a1] = 0.5 * model_concentrations[0][i][idx_A];
    state.variables_[i][_a2] = 0.5 * model_concentrations[0][i][idx_A];
    state.variables_[i][_b] = model_concentrations[0][i][idx_B];
    state.variables_[i][_c] = model_concentrations[0][i][idx_C];
    state.conditions_[i].temperature_ = temperatures[i];
    state.conditions_[i].pressure_ = pressures[i];
    state.conditions_[i].air_density_ = calculate_air_density_mol_m3(pressures[i], temperatures[i]);
  }

  std::vector<double> times;
  times.push_back(0);
  for (size_t i_time = 1; i_time < nsteps; ++i_time)
  {
    solver.CalculateRateConstants(state);
    prepare_for_solve(state);
    // Model results
    auto result = solver.Solve(time_step, state);
    postpare_for_solve(state);
    EXPECT_EQ(result.state_, (micm::SolverState::Converged));
    for (std::size_t i = 0; i < NUM_CELLS; ++i)
    {
      model_concentrations[i_time][i][idx_A] = state.variables_[i][_a1] + state.variables_[i][_a2];
      model_concentrations[i_time][i][idx_B] = state.variables_[i][_b];
      model_concentrations[i_time][i][idx_C] = state.variables_[i][_c];
    }

    // Analytical results
    double time = i_time * time_step;
    times.push_back(time);

    for (std::size_t i = 0; i < NUM_CELLS; ++i)
    {
      double initial_A = analytical_concentrations[0][i][idx_A];
      analytical_concentrations[i_time][i][idx_A] = initial_A * std::exp(-k1[i] * time);
      analytical_concentrations[i_time][i][idx_B] =
          initial_A * (k1[i] / (k2[i] - k1[i])) * (std::exp(-k1[i] * time) - std::exp(-k2[i] * time));
      analytical_concentrations[i_time][i][idx_C] =
          initial_A * (1.0 + (k1[i] * std::exp(-k2[i] * time) - k2[i] * std::exp(-k1[i] * time)) / (k2[i] - k1[i]));
    }
  }

  std::vector<std::string> header = { "time", "cell", "A", "B", "C" };
  writeCSV2D(test_label + "_stiff_model_concentrations.csv", header, model_concentrations, times);
  writeCSV2D(test_label + "_stiff_analytical_concentrations.csv", header, analytical_concentrations, times);

  for (std::size_t i_time = 1; i_time < model_concentrations.size(); ++i_time)
  {
    for (std::size_t i_cell = 0; i_cell < model_concentrations[i_time].size(); ++i_cell)
    {
      EXPECT_NEAR(
          model_concentrations[i_time][i_cell][idx_A], analytical_concentrations[i_time][i_cell][idx_A], absolute_tolerances)
          << "Arrays differ at index (" << i_time << ", " << i_cell << ", " << 0 << ") for " << test_label;
      EXPECT_NEAR(
          model_concentrations[i_time][i_cell][idx_B], analytical_concentrations[i_time][i_cell][idx_B], absolute_tolerances)
          << "Arrays differ at index (" << i_time << ", " << i_cell << ", " << 1 << ") for " << test_label;
      EXPECT_NEAR(
          model_concentrations[i_time][i_cell][idx_C], analytical_concentrations[i_time][i_cell][idx_C], absolute_tolerances)
          << "Arrays differ at index (" << i_time << ", " << i_cell << ", " << 2 << ") for " << test_label;
    }
  }
}

///////////////////////////////
// Specific analytical tests //
///////////////////////////////

template<class BuilderPolicy>
void test_analytical_troe(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-10,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TroeRateConstant({ .k0_A_ = 4.0e-11 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::TroeRateConstant({ .k0_A_ = 1.2e-3,
                                                                   .k0_B_ = 1.6,
                                                                   .k0_C_ = 3,
                                                                   .kinf_A_ = 136,
                                                                   .kinf_B_ = 5,
                                                                   .kinf_C_ = 24,
                                                                   .Fc_ = 0.9,
                                                                   .N_ = 0.8 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_system<BuilderPolicy>(
      "troe",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k_0 = 4.0e-11;
        double k_inf = 1;
        return k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.6, 1.0 / (1.0 + std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k_0 = 1.2e-3 * std::exp(3.0 / temperature) * std::pow(temperature / 300.0, 1.6);
        double k_inf = 136.0 * std::exp(24.0 / temperature) * std::pow(temperature / 300.0, 5.0);
        return k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.9, 0.8 / (0.8 + std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_stiff_troe(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-5,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TroeRateConstant({ .k0_A_ = 4.0e-11 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TroeRateConstant({ .k0_A_ = 4.0e-11 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::TroeRateConstant({ .k0_A_ = 1.2e-3,
                                                                   .k0_B_ = 1.6,
                                                                   .k0_C_ = 3,
                                                                   .kinf_A_ = 136,
                                                                   .kinf_B_ = 5,
                                                                   .kinf_C_ = 24,
                                                                   .Fc_ = 0.9,
                                                                   .N_ = 0.8 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(a2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(a1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_stiff_system<BuilderPolicy>(
      "troe",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k_0 = 4.0e-11;
        double k_inf = 1;
        return k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.6, 1.0 / (1.0 + std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k_0 = 1.2e-3 * std::exp(3.0 / temperature) * std::pow(temperature / 300.0, 1.6);
        double k_inf = 136.0 * std::exp(24.0 / temperature) * std::pow(temperature / 300.0, 5.0);
        return k_0 * air_density / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_photolysis(
    BuilderPolicy builder,
    double absolute_tolerances = 2e-6,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "photoA" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "photoB" }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  std::unordered_map<std::string, std::vector<double>> custom_parameters = {
    { "photoA", std::vector<double>(NUM_CELLS, 2e-3) }, { "photoB", std::vector<double>(NUM_CELLS, 3e-3) }
  };

  test_simple_system<BuilderPolicy>(
      "photolysis",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        return 2e-3;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        return 3e-3;
      },
      prepare_for_solve,
      postpare_for_solve,
      custom_parameters);
}

template<class BuilderPolicy>
void test_analytical_stiff_photolysis(
    BuilderPolicy builder,
    double absolute_tolerances = 2e-5,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "photoA1B" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "photoA2B" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "photoB" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(a2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(a1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  std::unordered_map<std::string, std::vector<double>> custom_parameters = {
    { "photoA1B", std::vector<double>(NUM_CELLS, 2e-3) },
    { "photoA2B", std::vector<double>(NUM_CELLS, 2e-3) },
    { "photoB", std::vector<double>(NUM_CELLS, 3e-3) }
  };

  test_simple_stiff_system<BuilderPolicy>(
      "photolysis",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        return 2e-3;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        return 3e-3;
      },
      prepare_for_solve,
      postpare_for_solve,
      custom_parameters);
}

template<class BuilderPolicy>
void test_analytical_ternary_chemical_activation(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-08,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 4.0e-5, .kinf_A_ = 1 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2e-3,
                                                                                        .k0_B_ = 1.6,
                                                                                        .k0_C_ = 3,
                                                                                        .kinf_A_ = 136,
                                                                                        .kinf_B_ = 5,
                                                                                        .kinf_C_ = 24,
                                                                                        .Fc_ = 0.9,
                                                                                        .N_ = 0.8 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_system<BuilderPolicy>(
      "ternary_chemical_activation",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k_0 = 4.0e-5;
        double k_inf = 1;
        return k_0 / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.6, 1.0 / (1.0 + (1.0 / 1.0) * std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k_0 = 1.2e-3 * std::exp(3.0 / temperature) * std::pow(temperature / 300.0, 1.6);
        double k_inf = 136.0 * std::exp(24.0 / temperature) * std::pow(temperature / 300.0, 5.0);
        return k_0 / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_stiff_ternary_chemical_activation(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-6,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 4.0e-5, .kinf_A_ = 1 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 4.0e-5, .kinf_A_ = 1 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2e-3,
                                                                                        .k0_B_ = 1.6,
                                                                                        .k0_C_ = 3,
                                                                                        .kinf_A_ = 136,
                                                                                        .kinf_B_ = 5,
                                                                                        .kinf_C_ = 24,
                                                                                        .Fc_ = 0.9,
                                                                                        .N_ = 0.8 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(a2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(a1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_stiff_system<BuilderPolicy>(
      "ternary_chemical_activation",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k_0 = 4.0e-5;
        double k_inf = 1;
        return k_0 / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.6, 1.0 / (1.0 + std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k_0 = 1.2e-3 * std::exp(3.0 / temperature) * std::pow(temperature / 300.0, 1.6);
        double k_inf = 136.0 * std::exp(24.0 / temperature) * std::pow(temperature / 300.0, 5.0);
        return k_0 / (1.0 + k_0 * air_density / k_inf) *
               std::pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * std::pow(std::log10(k_0 * air_density / k_inf), 2)));
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_tunneling(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-8,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TunnelingRateConstant({ .A_ = 4.0e-3 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::TunnelingRateConstant({ .A_ = 1.2e-4, .B_ = 167, .C_ = 1.0e8 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_system<BuilderPolicy>(
      "tunneling",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k1 = 4.0e-3;
        return k1;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k2 = 1.2e-4 * std::exp(-167 / temperature + 1.0e8 / std::pow(temperature, 3));
        return k2;
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_stiff_tunneling(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-6,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TunnelingRateConstant({ .A_ = 4.0e-3 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::TunnelingRateConstant({ .A_ = 4.0e-3 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::TunnelingRateConstant({ .A_ = 1.2e-4, .B_ = 167, .C_ = 1.0e8 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(a2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(a1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_stiff_system<BuilderPolicy>(
      "tunneling",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k1 = 4.0e-3;
        return k1;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k2 = 1.2e-4 * std::exp(-167 / temperature + 1.0e8 / std::pow(temperature, 3));
        return k2;
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_arrhenius(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-9,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e-3, .C_ = 50 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 =
      micm::ChemicalReactionBuilder()
          .SetReactants({ b })
          .SetProducts({ micm::Yield(c, 1) })
          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 1.2e-4, .B_ = 7, .C_ = 75, .D_ = 50, .E_ = 0.5 }))
          .SetPhaseName("gas")
          .Build();

  auto processes = std::vector<micm::Process>{ r1, r2 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_system<BuilderPolicy>(
      "arrhenius",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k1 = 4.0e-3 * std::exp(50 / temperature);
        return k1;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k2 = 1.2e-4 * std::exp(75 / temperature) * std::pow(temperature / 50, 7) * (1.0 + 0.5 * pressure);
        return k2;
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_stiff_arrhenius(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-6,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e-3, .C_ = 50 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e-3, .C_ = 50 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 =
      micm::ChemicalReactionBuilder()
          .SetReactants({ b })
          .SetProducts({ micm::Yield(c, 1) })
          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 1.2e-4, .B_ = 1.6, .C_ = 75, .D_ = 50, .E_ = 0.5 }))
          .SetPhaseName("gas")
          .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(a2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(a1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_stiff_system<BuilderPolicy>(
      "arrhenius",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double k1 = 4.0e-3 * std::exp(50 / temperature);
        return k1;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double k2 = 1.2e-4 * std::exp(75 / temperature) * std::pow(temperature / 50, 1.6) * (1.0 + 0.5 * pressure);
        return k2;
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_branched(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-13,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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
      micm::ChemicalReactionBuilder()
          .SetReactants({ a })
          .SetProducts({ micm::Yield(b, 1) })
          .SetRateConstant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy,
                                                        .X_ = 1e-4,
                                                        .Y_ = 204.3,
                                                        .a0_ = 1.0e-3,
                                                        .n_ = 2 }))
          .SetPhaseName("gas")
          .Build();

  micm::Process r2 =
      micm::ChemicalReactionBuilder()
          .SetReactants({ b })
          .SetProducts({ micm::Yield(c, 1) })
          .SetRateConstant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Nitrate,
                                                        .X_ = 1e-4,
                                                        .Y_ = 204.3,
                                                        .a0_ = 1.0e-3,
                                                        .n_ = 2 }))
          .SetPhaseName("gas")
          .Build();

  auto processes = std::vector<micm::Process>{ r1, r2 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_system<BuilderPolicy>(
      "branched",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double air_dens_n_cm3 = air_density * micm::constants::AVOGADRO_CONSTANT * 1.0e-6;
        double a_ = 2.0e-22 * std::exp(2) * 2.45e19;
        double b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
        double z =
            a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
        a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
        b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
        double A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

        double k1 = 1e-4 * std::exp(-204.3 / temperature) * (z / (z + A));
        return k1;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double air_dens_n_cm3 = air_density * micm::constants::AVOGADRO_CONSTANT * 1.0e-6;
        double a_ = 2.0e-22 * std::exp(2) * 2.45e19;
        double b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
        double z =
            a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
        a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
        b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
        double A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

        double k2 = 1e-4 * std::exp(-204.3 / temperature) * (A / (z + A));
        return k2;
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_stiff_branched(
    BuilderPolicy builder,
    double absolute_tolerances = 1e-6,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
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
      micm::ChemicalReactionBuilder()
          .SetReactants({ a1 })
          .SetProducts({ micm::Yield(b, 1) })
          .SetRateConstant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy,
                                                        .X_ = 1e-4,
                                                        .Y_ = 204.3,
                                                        .a0_ = 1.0e-3,
                                                        .n_ = 2 }))
          .SetPhaseName("gas")
          .Build();

  micm::Process r2 =
      micm::ChemicalReactionBuilder()
          .SetReactants({ a2 })
          .SetProducts({ micm::Yield(b, 1) })
          .SetRateConstant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy,
                                                        .X_ = 1e-4,
                                                        .Y_ = 204.3,
                                                        .a0_ = 1.0e-3,
                                                        .n_ = 2 }))
          .SetPhaseName("gas")
          .Build();

  micm::Process r3 =
      micm::ChemicalReactionBuilder()
          .SetReactants({ b })
          .SetProducts({ micm::Yield(c, 1) })
          .SetRateConstant(micm::BranchedRateConstant({ .branch_ = micm::BranchedRateConstantParameters::Branch::Nitrate,
                                                        .X_ = 1e-4,
                                                        .Y_ = 204.3,
                                                        .a0_ = 1.0e-3,
                                                        .n_ = 2 }))
          .SetPhaseName("gas")
          .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(a2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2 })
                         .SetProducts({ micm::Yield(a1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.9 * 4.0e10 }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5 };
  builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase })).SetReactions(processes);

  test_simple_stiff_system<BuilderPolicy>(
      "branched",
      builder,
      absolute_tolerances,
      [](double temperature, double pressure, double air_density)
      {
        // A->B reaction rate
        double air_dens_n_cm3 = air_density * micm::constants::AVOGADRO_CONSTANT * 1.0e-6;
        double a_ = 2.0e-22 * std::exp(2) * 2.45e19;
        double b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
        double z =
            a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
        a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
        b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
        double A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

        double k1 = 1e-4 * std::exp(-204.3 / temperature) * (z / (z + A));
        return k1;
      },
      [](double temperature, double pressure, double air_density)
      {
        // B->C reaction rate
        double air_dens_n_cm3 = air_density * micm::constants::AVOGADRO_CONSTANT * 1.0e-6;
        double a_ = 2.0e-22 * std::exp(2) * 2.45e19;
        double b_ = 0.43 * std::pow((293.0 / 298.0), -8.0);
        double z =
            a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
        a_ = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
        b_ = 0.43 * std::pow((temperature / 298.0), -8.0);
        double A = a_ / (1.0 + a_ / b_) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_ / b_), 2)));

        double k2 = 1e-4 * std::exp(-204.3 / temperature) * (A / (z + A));
        return k2;
      },
      prepare_for_solve,
      postpare_for_solve);
}

template<class BuilderPolicy>
void test_analytical_robertson(
    BuilderPolicy builder,
    double relative_tolerance = 1e-8,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
{
  /*
   * A -> B, k1 = 0.04
   * B + B -> C + B, k2 = 3e7
   * B + C -> A + C, k3 = 1e4
   *
   * this problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b, b })
                         .SetProducts({ micm::Yield(b, 1), micm::Yield(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b, c })
                         .SetProducts({ micm::Yield(a, 1), micm::Yield(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3 };
  auto solver = builder.SetReorderState(false)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions(processes)
                    .Build();

  double temperature = 272.5;
  double pressure = 101253.3;
  double air_density = 1e6;

  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1e-10);
  state.SetAbsoluteTolerances(std::vector<double>(3, state.relative_tolerance_ * 1e-2));

  double k1 = 0.04;
  double k2 = 3e7;
  double k3 = 1e4;

  state.SetCustomRateParameter("r1", k1);
  state.SetCustomRateParameter("r2", k2);
  state.SetCustomRateParameter("r3", k3);

  constexpr size_t N = 12;

  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(3));
  std::vector<std::vector<double>> analytical_concentrations(N + 1, std::vector<double>(3));

  model_concentrations[0] = { 1, 0, 0 };

  analytical_concentrations = { { 1, 0, 0 },
                                { 0.9664597373330035E+00, 0.3074626578578675E-04, 0.3350951640121071E-01 },
                                { 0.8413699238414729E+00, 0.1623390937990473E-04, 0.1586138422491472E+00 },
                                { 0.6172348823960878E+00, 0.6153591274639123E-05, 0.3827589640126376E+00 },
                                { 0.3368745306607069E+00, 0.2013702318261393E-05, 0.6631234556369748E+00 },
                                { 0.1073004285378040E+00, 0.4800166972571660E-06, 0.8926990914454987E+00 },
                                { 0.1786592114209946E-01, 0.7274751468436319E-07, 0.9821340061103859E+00 },
                                { 0.2031483924973415E-02, 0.8142277783356159E-08, 0.9979685079327488E+00 },
                                { 0.2076093439016395E-03, 0.8306077485067610E-09, 0.9997923898254906E+00 },
                                { 0.2082417512179460E-04, 0.8329841429908955E-10, 0.9999791757415798E+00 },
                                { 0.2083229471647004E-05, 0.8332935037760723E-11, 0.9999979167621954E+00 },
                                { 0.2083328471883087E-06, 0.8333315602809495E-12, 0.9999997916663195E+00 },
                                { 0.2083340149701284E-07, 0.8333360770334744E-13, 0.9999999791665152E+00 } };

  state.variables_[0] = model_concentrations[0];
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  double time_step = 1.0;
  std::vector<double> times;
  times.push_back(0);
  solver.CalculateRateConstants(state);
  for (size_t i_time = 0; i_time < N; ++i_time)
  {
    double solve_time = time_step + i_time * time_step;
    times.push_back(solve_time);
    prepare_for_solve(state);
    // Model results
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      auto result = solver.Solve(time_step - actual_solve, state);
      actual_solve += result.final_time_;
    }
    postpare_for_solve(state);
    model_concentrations[i_time + 1] = state.variables_[0];
    time_step *= 10;
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("robertson_model_concentrations.csv", header, model_concentrations, times);
  writeCSV("robertson_analytical_concentrations.csv", header, analytical_concentrations, times);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  // average of the starting concentration *1e-3;
  double absolute_tolerance = 0.3e-3;
  for (size_t i = 1; i < model_concentrations.size(); ++i)
  {
    double rel_error = relative_error(model_concentrations[i][_a], analytical_concentrations[i][0]);
    double abs_error = std::abs(model_concentrations[i][_a] - analytical_concentrations[i][0]);
    EXPECT_TRUE(abs_error < absolute_tolerance || rel_error < relative_tolerance)
        << "Arrays differ at index (" << i << ", " << 0 << ") with relative error " << rel_error << " and absolute error "
        << abs_error;

    rel_error = relative_error(model_concentrations[i][_b], analytical_concentrations[i][1]);
    abs_error = std::abs(model_concentrations[i][_b] - analytical_concentrations[i][1]);
    EXPECT_TRUE(abs_error < absolute_tolerance || rel_error < relative_tolerance)
        << "Arrays differ at index (" << i << ", " << 1 << ") with relative error " << rel_error << " and absolute error "
        << abs_error;

    rel_error = relative_error(model_concentrations[i][_c], analytical_concentrations[i][2]);
    abs_error = std::abs(model_concentrations[i][_c] - analytical_concentrations[i][2]);
    EXPECT_TRUE(abs_error < absolute_tolerance || rel_error < relative_tolerance)
        << "Arrays differ at index (" << i << ", " << 2 << ") with relative error " << rel_error << " and absolute error "
        << abs_error;
  }
}

template<class BuilderPolicy>
void test_analytical_oregonator(
    BuilderPolicy builder,
    double absolute_tolerance = 1e-8,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
{
  /*
   * This problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 144. It actually comes from Field and Noyes (1974)
   * A driver for a version of this is from here, but this needs a custom forcing and jacobian and so we tried to translate
   * it https://www.unige.ch/~hairer/testset/testset.html
   *
   * Field, R.J., Noyes, R.M., 1974. Oscillations in chemical systems. IV. Limit cycle behavior in a model of a real chemical
   * reaction. The Journal of Chemical Physics 60, 1877–1884. https://doi.org/10.1063/1.1681288
   *
   * In this paper, they give five equations that use X, Y, Z, A, P, Q, and B. P and Q are only produced and don't react.
   * They set A = B = [BrO3-] = 0.06. Those equations come to this
   *
   * Y -> X,      k1 = 1.34 * 0.06
   * X + Y -> P,  k2 = 1.6e9
   * X -> Z + 2X, k3 = 8e3*0.06
   * 2X -> Q,     k4 = 4e7
   * Z -> Y,      k5 = 1
   *
   * Through some other more complicatd math they simplifed to only 3 variables, X, Y, and Z, but I couldn't figure out how
   * to represent those equations
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   *
   * I don't understand the transfomrations. Multiplying the timestep by tau, and the concnetrations by the constants
   * in the paper give very similar values.
   */

  auto X = micm::Species("X");
  auto Y = micm::Species("Y");
  auto Z = micm::Species("Z");
  auto P = micm::Species("P");
  auto Q = micm::Species("Q");

  micm::Phase gas_phase{ std::vector<micm::Species>{ X, Y, Z, P, Q } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ Y })
                         .SetProducts({ micm::Yield(X, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ X, Y })
                         .SetProducts({ micm::Yield(P, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ X })
                         .SetProducts({ micm::Yield(Z, 1), micm::Yield(X, 2) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ X, X })
                         .SetProducts({ micm::Yield(Q, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r4" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ Z })
                         .SetProducts({ micm::Yield(Y, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r5" }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5 };
  auto solver = builder.SetReorderState(false)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions(processes)
                    .Build();

  double tau = 0.1610;
  double time_step = 30 * tau;
  size_t N = 12;

  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(5));
  std::vector<std::vector<double>> analytical_concentrations(N + 1, std::vector<double>(5));

  // ignore P and Q, the last two zeros
  // the initial concentrations given at https://www.unige.ch/~hairer/testset/testset.html
  // are for the dimensionless variables, so we need to convert them to the actual concentrations
  // X = alpha
  // Y = eta
  // Z = rho
  double alpha_const = 5.025e-11;
  double eta_const = 3e-7;
  double rho_const = 2.412e-8;
  model_concentrations[0] = { 1 * alpha_const, 2 * eta_const, 3 * rho_const, 0, 0 };

  // ignore P and Q, the last two zeros
  analytical_concentrations = {
    { 1, 2, 3, 0, 0 },
    { 0.1000661467180497E+01, 0.1512778937348249E+04, 0.1035854312767229E+05, 0, 0 },
    { 0.1000874625199626E+01, 0.1144336972384497E+04, 0.8372149966624639E+02, 0, 0 },
    { 0.1001890368438751E+01, 0.5299926232295553E+03, 0.1662279579042420E+01, 0, 0 },
    { 0.1004118022612645E+01, 0.2438326079910346E+03, 0.1008822224048647E+01, 0, 0 },
    { 0.1008995416634061E+01, 0.1121664388662539E+03, 0.1007783229065319E+01, 0, 0 },
    { 0.1019763472537298E+01, 0.5159761322947535E+02, 0.1016985778956374E+01, 0, 0 },
    { 0.1043985088527474E+01, 0.2373442027531524E+02, 0.1037691843544522E+01, 0, 0 },
    { 0.1100849071667922E+01, 0.1091533805469020E+02, 0.1085831969810860E+01, 0, 0 },
    { 0.1249102130020572E+01, 0.5013945178605446E+01, 0.1208326626237875E+01, 0, 0 },
    { 0.1779724751937019E+01, 0.2281852385542403E+01, 0.1613754023671725E+01, 0, 0 },
    { 0.1000889326903503E+01, 0.1125438585746596E+04, 0.1641049483777168E+05, 0, 0 },
    { 0.1000814870318523E+01, 0.1228178521549889E+04, 0.1320554942846513E+03, 0, 0 },
  };

  for (auto& row : analytical_concentrations)
  {
    row[0] *= alpha_const;
    row[1] *= eta_const;
    row[2] *= rho_const;
  }

  auto state = solver.GetState(1);

  state.SetRelativeTolerance(1e-6);
  state.SetAbsoluteTolerances(std::vector<double>(5, state.relative_tolerance_ * 1e-6));

  state.SetCustomRateParameter("r1", 1.34 * 0.06);
  state.SetCustomRateParameter("r2", 1.6e9);
  state.SetCustomRateParameter("r3", 8e3 * 0.06);
  state.SetCustomRateParameter("r4", 4e7);
  state.SetCustomRateParameter("r5", 1);

  state.variables_[0] = model_concentrations[0];
  solver.CalculateRateConstants(state);
  prepare_for_solve(state);

  std::vector<double> times;
  times.push_back(0);
  for (size_t i_time = 0; i_time < N; ++i_time)
  {
    double solve_time = time_step + i_time * time_step;
    times.push_back(solve_time);
    // Model results
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      auto result = solver.Solve(time_step - actual_solve, state);
      actual_solve += result.final_time_;
    }
    postpare_for_solve(state);
    model_concentrations[i_time + 1] = state.variables_[0];
  }

  std::vector<std::string> header = { "time", "X", "Y", "Z", "P", "Q" };
  writeCSV("oregonator_model_concentrations.csv", header, model_concentrations, times);
  std::vector<double> an_times;
  an_times.push_back(0);
  for (int i = 1; i <= 12; ++i)
  {
    an_times.push_back(time_step * i);
  }
  writeCSV("oregonator_analytical_concentrations.csv", header, analytical_concentrations, an_times);

  auto map = state.variable_map_;

  size_t _x = map.at("X");
  size_t _y = map.at("Y");
  size_t _z = map.at("Z");

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_x], analytical_concentrations[i][0], absolute_tolerance);
    EXPECT_NEAR(model_concentrations[i][_y], analytical_concentrations[i][1], absolute_tolerance);
    EXPECT_NEAR(model_concentrations[i][_z], analytical_concentrations[i][2], absolute_tolerance);
  }
}

template<class BuilderPolicy>
void test_analytical_hires(
    BuilderPolicy builder,
    double absolute_tolerance = 1e-8,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
{
  /*
   * This problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 144
   *
   * From the forcing function these equations were made
   * y0 -> y1,                                  k1 = 1.71
   * y1 -> (0.43 / 8.75 )y0 + (8.32 / 8.75)y3,  k2 = 8.75
   * y2 -> (8.32 /10.03)y0 + (1.71/10.03)y3,    k3 = 10.03
   *    -> y0,                                  k4 = 0.0007
   * y3 -> (0.43/1.12)y2 + (0.69/1.12)y5,       k5 = 1.12
   * y4 -> (0.035/1.745)y2 + (1.71/1.745)y5,    k6 = 1.745
   * y5 -> y4                                   k7 = 0.43
   * y6 -> (0.43/1.81)y4 + (0.69/1.81)y5 + y7   k8 = 1.81
   * y5 + y7 -> y6                              k9 = 280.0
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto y0 = micm::Species("Y0");
  auto y1 = micm::Species("Y1");
  auto y2 = micm::Species("Y2");
  auto y3 = micm::Species("Y3");
  auto y4 = micm::Species("Y4");
  auto y5 = micm::Species("Y5");
  auto y6 = micm::Species("Y6");
  auto y7 = micm::Species("Y7");

  micm::Phase gas_phase{ std::vector<micm::Species>{ y0, y1, y2, y3, y4, y5, y6, y7 } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y0 })
                         .SetProducts({ micm::Yield(y1, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y1 })
                         .SetProducts({ micm::Yield(y0, 0.43 / 8.75), micm::Yield(y3, 8.32 / 8.75) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y2 })
                         .SetProducts({ micm::Yield(y0, 8.32 / 10.03), micm::Yield(y3, 1.71 / 10.03) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetProducts({ micm::Yield(y0, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r4" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r5 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y3 })
                         .SetProducts({ micm::Yield(y2, 0.43 / 1.12), micm::Yield(y5, 0.69 / 1.12) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r5" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r6 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y4 })
                         .SetProducts({ micm::Yield(y2, 0.035 / 1.745), micm::Yield(y5, 1.71 / 1.745) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r6" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r7 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y5 })
                         .SetProducts({ micm::Yield(y4, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r7" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r8 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y6 })
                         .SetProducts({ micm::Yield(y4, 0.43 / 1.81), micm::Yield(y5, 0.69 / 1.81), micm::Yield(y7, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r8" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r9 = micm::ChemicalReactionBuilder()
                         .SetReactants({ y5, y7 })
                         .SetProducts({ micm::Yield(y6, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r9" }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6, r7, r8, r9 };
  auto solver = builder.SetReorderState(false)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions(processes)
                    .Build();

  size_t N = 2;
  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(8));
  std::vector<std::vector<double>> analytical_concentrations(N + 1, std::vector<double>(8));

  model_concentrations[0] = { 1, 0, 0, 0, 0, 0, 0, 0.0057 };

  analytical_concentrations = {
    { 1, 0, 0, 0, 0, 0, 0, 0.0057 },
    { 0.000737131257332567,
      0.000144248572631618,
      0.000058887297409676,
      0.001175651343283149,
      0.002386356198831330,
      0.006238968252742796,
      0.002849998395185769,
      0.002850001604814231 },
    { 0.000670305503581864,
      0.000130996846986347,
      0.000046862231597733,
      0.001044668020551705,
      0.000594883830951485,
      0.001399628833942774,
      0.001014492757718480,
      0.004685507242281520 },
  };

  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1e-6);
  state.SetAbsoluteTolerances(std::vector<double>(8, state.relative_tolerance_ * 1e-2));

  state.SetCustomRateParameter("r1", 1.71);
  state.SetCustomRateParameter("r2", 8.75);
  state.SetCustomRateParameter("r3", 10.03);
  state.SetCustomRateParameter("r4", 0.0007);
  state.SetCustomRateParameter("r5", 1.12);
  state.SetCustomRateParameter("r6", 1.745);
  state.SetCustomRateParameter("r7", 0.43);
  state.SetCustomRateParameter("r8", 1.81);
  state.SetCustomRateParameter("r9", 280.0);

  state.variables_[0] = model_concentrations[0];
  solver.CalculateRateConstants(state);
  prepare_for_solve(state);

  std::vector<double> times;
  times.push_back(0);
  double time_step = 321.8122;
  for (size_t i_time = 0; i_time < N; ++i_time)
  {
    double solve_time = time_step + i_time * time_step;
    times.push_back(solve_time);
    // Model results
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      auto result = solver.Solve(time_step - actual_solve, state);
      actual_solve += result.final_time_;
    }
    postpare_for_solve(state);
    model_concentrations[i_time + 1] = state.variables_[0];
    time_step += 100;
  }

  std::vector<std::string> header = { "time", "y0", "y1", "y2", "y3", "y4", "y5", "y6", "y7" };
  writeCSV("hires_model_concentrations.csv", header, model_concentrations, times);
  writeCSV("hires_analytical_concentrations.csv", header, analytical_concentrations, times);

  auto map = state.variable_map_;

  size_t _y0 = map.at("Y0");
  size_t _y1 = map.at("Y1");
  size_t _y2 = map.at("Y2");
  size_t _y3 = map.at("Y3");
  size_t _y4 = map.at("Y4");
  size_t _y5 = map.at("Y5");
  size_t _y6 = map.at("Y6");
  size_t _y7 = map.at("Y7");

  for (size_t i = 1; i < model_concentrations.size(); ++i)
  {
    EXPECT_NEAR(model_concentrations[i][_y0], analytical_concentrations[i][0], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 0 << ")";
    EXPECT_NEAR(model_concentrations[i][_y1], analytical_concentrations[i][1], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 1 << ")";
    EXPECT_NEAR(model_concentrations[i][_y2], analytical_concentrations[i][2], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 2 << ")";
    EXPECT_NEAR(model_concentrations[i][_y3], analytical_concentrations[i][3], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 3 << ")";
    EXPECT_NEAR(model_concentrations[i][_y4], analytical_concentrations[i][4], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 4 << ")";
    EXPECT_NEAR(model_concentrations[i][_y5], analytical_concentrations[i][5], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 5 << ")";
    EXPECT_NEAR(model_concentrations[i][_y6], analytical_concentrations[i][6], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 6 << ")";
    EXPECT_NEAR(model_concentrations[i][_y7], analytical_concentrations[i][7], absolute_tolerance)
        << "Arrays differ at index (" << i << ", " << 7 << ")";
  }
}

template<class BuilderPolicy>
void test_analytical_e5(
    BuilderPolicy builder,
    double relative_tolerance = 1e-8,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {},
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) {})
{
  /*
   * A1 -> A2 + A3,  k1 = 7.89e-10
   * A2 + A3 -> A5,  k2 = 1.13e9
   * A1 + A3 -> A4,  k3 = 1.1e7
   * A4 -> A3 + A6,  k4 = 1.13e3
   *
   * this problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * full equations retrieved from here: https://archimede.uniba.it/~testset/report/e5.pdf
   *
   * originally described here
   *
   * Enright, W.H., Hull, T.E., Lindberg, B., 1975. Comparing numerical methods for stiff systems of O.D.E:s. BIT 15, 10–48.
   * https://doi.org/10.1007/BF01932994
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto a1 = micm::Species("A1");
  auto a2 = micm::Species("A2");
  auto a3 = micm::Species("A3");
  auto a4 = micm::Species("A4");
  auto a5 = micm::Species("A5");
  auto a6 = micm::Species("A6");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a1, a2, a3, a4, a5, a6 } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1 })
                         .SetProducts({ micm::Yield(a2, 1), micm::Yield(a3, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a2, a3 })
                         .SetProducts({ micm::Yield(a5, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a1, a3 })
                         .SetProducts({ micm::Yield(a4, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a4 })
                         .SetProducts({ micm::Yield(a3, 1), micm::Yield(a6, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r4" }))
                         .SetPhaseName("gas")
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3, r4 };
  auto solver = builder.SetReorderState(false)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions(processes)
                    .Build();

  size_t N = 7;

  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(6));
  std::vector<std::vector<double>> analytical_concentrations(N + 1, std::vector<double>(6));

  model_concentrations[0] = { 1.76e-3, 0, 0, 0, 0, 0 };

  // ignore the concentration of A5 and A6
  analytical_concentrations = {
    { 1.76e-3, 0, 0, 0, 0, 0 },
    { 1.7599259497677897058e-003, 1.3846281519376516449e-011, 7.6370038530073911180e-013, 1.3082581134075777338e-011, 0, 0 },
    { 1.6180769999072942552e-003, 1.3822370304983735443e-010, 8.2515735006838336088e-012, 1.2997212954915352082e-010, 0, 0 },
    { 7.4813208224292220114e-006, 2.3734781561205975019e-012, 2.2123586689581663654e-012, 1.6111948716243113653e-013, 0, 0 },
    { 4.7150333630401632232e-010, 1.8188895860807021729e-014, 1.8188812376786725407e-014, 8.3484020296321693074e-020, 0, 0 },
    { 3.1317148329356996037e-014, 1.4840957952870064294e-016, 1.4840957948345691466e-016, 4.5243728279782625194e-026, 0, 0 },
    { 3.8139035189787091771e-049, 1.0192582567660293322e-020, 1.0192582567660293322e-020, 3.7844935507486221171e-065, 0, 0 },
    { 0.0000000000000000000e-000, 8.8612334976263783420e-023, 8.8612334976263783421e-023, 0.0000000000000000000e-000, 0, 0 }
  };

  auto state = solver.GetState(1);

  state.SetRelativeTolerance(1e-13);
  state.SetAbsoluteTolerances(std::vector<double>(6, 1e-17));
  auto atol = state.absolute_tolerance_;
  atol[0] = 1e-7;
  atol[4] = 1e-7;
  atol[5] = 1e-7;
  state.SetAbsoluteTolerances(atol);

  state.SetCustomRateParameter("r1", 7.89e-10);
  state.SetCustomRateParameter("r2", 1.13e9);
  state.SetCustomRateParameter("r3", 1.1e7);
  state.SetCustomRateParameter("r4", 1.13e3);

  state.variables_[0] = model_concentrations[0];
  solver.CalculateRateConstants(state);
  prepare_for_solve(state);

  std::vector<double> times;
  times.push_back(0);
  double time_step = 10;
  for (size_t i_time = 0; i_time < N; ++i_time)
  {
    double solve_time = time_step + i_time * time_step;
    times.push_back(solve_time);
    // Model results
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      auto result = solver.Solve(time_step - actual_solve, state);
      actual_solve += result.final_time_;
    }
    postpare_for_solve(state);
    model_concentrations[i_time + 1] = state.variables_[0];
    time_step *= 100;
  }

  std::vector<std::string> header = { "time", "a1", "a2", "a3", "a4", "a5", "a6" };
  writeCSV("e5_model_concentrations.csv", header, model_concentrations, times);
  writeCSV("e5_analytical_concentrations.csv", header, analytical_concentrations, times);

  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    // ignore the concentration of A5 and A6
    double absolute_tolerance = 1e-6;
    double rel_error = relative_error(model_concentrations[i][0], analytical_concentrations[i][0]);
    double abs_error = std::abs(model_concentrations[i][0] - analytical_concentrations[i][0]);
    EXPECT_TRUE(abs_error < absolute_tolerance || rel_error < relative_tolerance)
        << "Arrays differ at index (" << i << ", " << 0 << ") with relative error " << rel_error << " and absolute error "
        << abs_error;

    absolute_tolerance = 1e-13;
    rel_error = relative_error(model_concentrations[i][1], analytical_concentrations[i][1]);
    abs_error = std::abs(model_concentrations[i][1] - analytical_concentrations[i][1]);
    EXPECT_TRUE(abs_error < absolute_tolerance || rel_error < relative_tolerance)
        << "Arrays differ at index (" << i << ", " << 1 << ") with relative error " << rel_error << " and absolute error "
        << abs_error;

    absolute_tolerance = 1e-13;
    rel_error = relative_error(model_concentrations[i][2], analytical_concentrations[i][2]);
    abs_error = std::abs(model_concentrations[i][2] - analytical_concentrations[i][2]);
    EXPECT_TRUE(abs_error < absolute_tolerance || rel_error < relative_tolerance)
        << "Arrays differ at index (" << i << ", " << 2 << ") with relative error " << rel_error << " and absolute error "
        << abs_error;

    absolute_tolerance = 1e-13;
    rel_error = relative_error(model_concentrations[i][3], analytical_concentrations[i][3]);
    abs_error = std::abs(model_concentrations[i][3] - analytical_concentrations[i][3]);
    EXPECT_TRUE(abs_error < absolute_tolerance || rel_error < relative_tolerance)
        << "Arrays differ at index (" << i << ", " << 3 << ") with relative error " << rel_error << " and absolute error "
        << abs_error;
  }
}
