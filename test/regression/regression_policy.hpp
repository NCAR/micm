#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

double calculate_air_density_mol_m3(double pressure, double temperature)
{
  return pressure / (micm::constants::GAS_CONSTANT * temperature);
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

std::pair<std::vector<std::string>, std::vector<std::vector<double>>> read_csv(const std::string& filename)
{
  std::ifstream file(filename);
  if (file.is_open())
  {
    std::vector<std::string> header;
    std::vector<std::vector<double>> data;

    // Read column headers
    std::string line;
    if (std::getline(file, line))
    {
      std::istringstream ss(line);
      std::string item;
      while (std::getline(ss, item, ','))
      {
        header.push_back(item);
      }
    }

    // Read data rows
    while (std::getline(file, line))
    {
      std::istringstream ss(line);
      std::string item;
      std::vector<double> row;
      while (std::getline(ss, item, ','))
      {
        row.push_back(std::stod(item));
      }

      data.push_back(row);
    }

    file.close();
    return { header, data };
  }
  else
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return {};
  }
}

template<class BuilderPolicy>
void test_flow_tube(
    BuilderPolicy builder,
    std::string expected_results_path,
    std::function<void(typename BuilderPolicy::StatePolicyType&)> prepare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) { },
    std::function<void(typename BuilderPolicy::StatePolicyType&)> postpare_for_solve =
        [](typename BuilderPolicy::StatePolicyType& state) { })
{
  /*
   * SOA1 ->                                , k1 = 0.01
   * SOA2 ->                                , k2 = 0.05
   * O3 + a-pinene -> 0.18 SOA1 + 0.09 SOA2 , k3 = Arrhenius(A=8.8e-17)
   *
   */

  auto soa1 = micm::Species("SOA1");
  auto soa2 = micm::Species("SOA2");
  auto apinene = micm::Species("a-pinene");
  auto o3 = micm::Species("O3");

  double MOLES_M3_TO_MOLECULES_CM3 = 1.0e-6 * 6.02214076e23;

  micm::Phase gas_phase{ std::vector<micm::Species>{ soa1, soa2, apinene, o3 } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ soa1 })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ soa2 })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ apinene, o3 })
                         .SetProducts({ micm::Yield(soa1, 0.18), micm::Yield(soa2, 0.09) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 8.8e-17 * MOLES_M3_TO_MOLECULES_CM3 }))
                         .SetPhase(gas_phase)
                         .Build();

  auto processes = std::vector<micm::Process>{ r1, r2, r3 };
  auto solver = builder.SetReorderState(false)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions(processes)
                    .Build();

  size_t N = 3600;

  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(4));

  auto state = solver.GetState();

  state.SetConcentration(soa1, 0);
  state.SetConcentration(soa2, 0);
  state.SetConcentration(apinene, 8e-8);
  state.SetConcentration(o3, 0.00002);

  state.SetCustomRateParameter("r1", 0.01);
  state.SetCustomRateParameter("r2", 0.05);

  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325;
  state.conditions_[0].air_density_ = calculate_air_density_mol_m3(101325, 298.15);

  solver.CalculateRateConstants(state);

  prepare_for_solve(state);

  model_concentrations[0] = state.variables_[0];

  std::vector<double> times;
  times.push_back(0);
  double time_step = 1;
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

  std::vector<std::string> header = { "time", "SOA1", "SOA2", "a-pinene", "O3" };
  // uncomment to write results to file if there are answer-breaking changes in a PR
  // writeCSV(expected_results_path, header, model_concentrations, times);

  // Check that the results we got match exactly to those we expect
  auto [header_out, data_out] = read_csv(expected_results_path);

  EXPECT_EQ(header, header_out);
  EXPECT_EQ(model_concentrations.size(), data_out.size());
  // data_out contains times in the first column and values in the rest
  // hence the minus one
  EXPECT_EQ(model_concentrations[0].size(), data_out[0].size() - 1);
  for (size_t i = 0; i < 10; ++i)
  {
    for (size_t j = 0; j < model_concentrations[i].size(); ++j)
    {
      EXPECT_NEAR(model_concentrations[i][j], data_out[i][j + 1], 1e-10)
          << "Arrays differ at index (" << i << ", " << j << ") with value " << model_concentrations[i][j]
          << " and expected value " << data_out[i][j + 1];
    }
  }
}