// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// GENERATED FILE. Do not edit. Re-run test/benchmark/import_ts1.py instead.
//
// The MOZART-TS1 mechanism (MZ327_TS1.2_20230307), as a benchmark mechanism.
// Read from the musica TS1 configuration, mechanism configuration version 1.0.0.
//
// 210 species and 547 reactions (Arrhenius 361, Photolysis 123, Surface 13, Troe 31, UserDefined 19).
// Use it to see how the linear solver and the rate calculations scale with
// mechanism size.
//
// No concentration in the CSV, set to zero: sink.
// No rate parameter in the CSV, set to zero: usr_CO_OH usr_DMS_OH.
#pragma once

#include <micm/CPU.hpp>
#include <micm/util/property_keys.hpp>

#include <algorithm>
#include <string_view>
#include <vector>

namespace bench
{
  namespace ts1
  {
    /// @brief Find a phase species by name. A surface reaction needs the phase
    ///        entry, because that is what carries the diffusion coefficient.
    inline const micm::PhaseSpecies& PhaseSpeciesByName(const micm::Phase& phase, std::string_view name)
    {
      for (const auto& phase_species : phase.phase_species_)
      {
        if (phase_species.species_.name_ == name)
        {
          return phase_species;
        }
      }
      throw std::runtime_error("ts1: no species named '" + std::string(name) + "' in the phase");
    }

    /// @brief Find a species by name in the phase. A reaction must use the
    ///        phase's own species object, because that is what carries the
    ///        third-body parameterization and the molecular weight.
    inline const micm::Species& SpeciesByName(const micm::Phase& phase, std::string_view name)
    {
      return PhaseSpeciesByName(phase, name).species_;
    }

    /// @brief The TS1 gas phase: 210 species.
    inline micm::Phase CreateGasPhase()
    {
      std::vector<micm::PhaseSpecies> phase_species;
      phase_species.reserve(210);

      {
        auto species = micm::Species("ALKNIT");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.133141 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BZOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.124135 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C6H5OOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.110109 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("COF2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("COFCL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HF");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("F");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BENZO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.159115 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BZOO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.123128 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("N2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("E90");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NH_5");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NH_50");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ST80_25");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PAN");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.121048 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MVK");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0700878 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MACROOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.120101 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SOAG0");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SOAG1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SOAG2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SOAG3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SOAG4");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa4_a1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa5_a1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa5_a2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa3_a1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa2_a1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa1_a1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa1_a2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa2_a2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa3_a2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("soa4_a2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ISOPNITB");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.147031980536874 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("SO3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("OCS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("S");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BCARYO2VBS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SF6");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("sink");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0010074 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MEK");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0721026 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MTERP");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("N2O5");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.108023495904642 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("HCN");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0270251 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SVOC");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ISOPNO3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.162118 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("RO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0890682 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PHENO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.175114 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BENZO2VBS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("IVOC");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERPNIT");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.215154758199771 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("ISOPO2VBS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("IVOCO2VBS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HNO3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0630123 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ACBZO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.137112 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3COOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0760498 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("SO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0640648 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MTERPO2VBS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH4");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0160406 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3CL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0504859 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3CO3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0750424 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C6H5O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.109102 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERPROD1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.168227 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HYAC");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0740762 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HPALD");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.116112 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TOLUO2VBS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H2O");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0180142 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.045803364407376784 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("EOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0780646 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NTERPOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.23118008751781 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("XYLEO2VBS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CCL4");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.153822 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CF2CLBR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.165365 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CF3BR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.14891 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CFC11");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.137368 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CFC113");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.187375 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CFC114");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.170921 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CFC115");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.154467 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CFC12");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.120913 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH2BR2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.173834 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3BR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0949372 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3CCL3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.133402 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0300061 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.079904 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BRCL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.115357 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BRO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0959034 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BRONO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.141909 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0354527 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CL2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0709054 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CL2O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.102904 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CLO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0514521 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CLONO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0974576 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HCOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0460246 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HBR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0809114 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HOBR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0969108 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HOCL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0524595 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("N");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0140067 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BIGENE");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0561032 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C2H4");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0280516 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C2H5O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0610578 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3COCHO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0720614 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3COCH3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0580768 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0159994 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("OCLO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0674515 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O1D");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0159994 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PHENO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.159115 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HCFC141B");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.116948 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HCFC142B");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.100494 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HCFC22");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0864679 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("DMS");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0621324 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C2H5OH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0460658 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HCL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0364601 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BEPOMUC");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.126109 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CHBR3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.25273 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H2402");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.259824 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0440098 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BZALD");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.106121 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BENZENE");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0781104 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C3H7O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0750836 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.047032 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BCARY");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.204343 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BIGALD");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0980982 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BIGALD2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0980982 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BIGALD3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0980982 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BIGALD4");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.112124 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BIGALK");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0721438 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H2O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0340136 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C2H5OOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0620652 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C2H6");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0300664 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C3H8");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0440922 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C3H6");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0420774 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH2O");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0300252 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3CN");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0410509 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C2H2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0260368 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3OH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.03204 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3OOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0480394 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CRESOL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.108136 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ENEO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.105109 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MACRO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.119093 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ISOPAO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.11712 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MALO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.115064 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ISOPBO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.11712 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MCO3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.101079 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MDIALO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.117079 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MEKO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.103094 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("EO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0770572 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("EO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0610578 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("GLYOXAL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0580064525191708 });
        // The configuration sets no diffusion limit for this species, so this is the
        // default the other surface species carry. It is not a measurement.
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("MPAN");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.147085 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NC4CH2OH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.147031980536874 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("ISOPOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.118127 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("GLYALD");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0600504 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.03307754409115893 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("HOCH2OO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0630314 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("H2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0020148 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HYDRALD");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.100113 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ISOP");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0681142 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NTERPO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.230232 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TOLO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.173141 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERP2O2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.199219 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYLENO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.187166 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERPO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.185234 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYLOLO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.203166 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PBZNIT");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.183118 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYLENES");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.106162 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.091083 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TOLUENE");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0921362 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.149119 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.150126 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERPROD2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.154201 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MEKOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.104101 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("MACR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0700878 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HONITR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.133362340623015 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("ISOPNITA");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.147031980536874 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("ISOPNOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.163125 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("IEPOX");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.118127 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ONITR");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.147031980536874 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("H2SO4");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0980784 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("N2O");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0440129 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NO3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.06186298085408304 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("OH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0170068 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PHENOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.176122 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("PHENOL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0941098 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYLOL");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.122161 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ROOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0900756 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("O3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0479982 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERPOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.186241 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TOLOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.174148 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYLENOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.188174 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("XYLOLOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.204173 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0280104 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3COOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0600504 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("CH3CHO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.044051 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BIGALD1");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0840724 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ALKOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.104143 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("DICARBO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.12909 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("BENZOOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.160122 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("ALKO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.103135 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("HO2NO2");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0790117 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("C3H7OOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.076091 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NH3");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.017031 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("TERP2OOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.200226 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("POOH");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0920904 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NOA");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.119074 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NC4CHO");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.14461174234895102 });
        phase_species.emplace_back(species, 1e-05);
      }
      {
        auto species = micm::Species("TEPOMUC");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.140134 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("NH4");
        species.SetProperty(micm::property_keys::MOLECULAR_WEIGHT, micm::Real{ 0.0 });
        phase_species.emplace_back(species);
      }
      {
        auto species = micm::Species("M");
        // Parameterized on air density, so it holds no state variable.
        species.SetThirdBody();
        phase_species.emplace_back(species);
      }

      return micm::Phase{ "gas", phase_species };
    }

    /// @brief The TS1 reactions: 547 processes.
    inline std::vector<micm::Process> CreateProcesses(const micm::Phase& gas_phase)
    {
      std::vector<micm::Process> processes;
      processes.reserve(547);

      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "TERPO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.95 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.025 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1204428.152, .B_ = 0.0, .C_ = 500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C3H7OOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "C3H7O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "N2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 873210.4102, .B_ = 0.0, .C_ = 220.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "CH2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 20475278.584, .B_ = 0.0, .C_ = -1600.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOLUO2VBS"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.13640049999999998 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0101005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.07630050000000001 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.2157005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.0738005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HCFC142B") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 782878.2988, .B_ = 0.0, .C_ = -1770.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CFC114"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 70459046.892, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4877934.0156, .B_ = 0.0, .C_ = 270.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HONITR"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ONITR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1204428.152, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPO2VBS"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0031005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0035005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0003005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0271005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.0474005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 127669.384112, .B_ = 0.0, .C_ = 1300.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3CO3"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1746420.8204, .B_ = 0.0, .C_ = 500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H2O"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 98160894.388, .B_ = 0.0, .C_ = 60.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 361328.4456, .B_ = 0.0, .C_ = 230.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "MTERP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MTERP"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MTERPO2VBS"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7226568.9120000005, .B_ = 0.0, .C_ = 440.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3BR"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 108398533.67999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "SVOC"), SpeciesByName(gas_phase, "OH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.5931004999999999 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.1534005 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.045900500000000004 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.008500500000000001 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.012800500000000001 } })
              .SetRateConstant(
                  micm::ArrheniusRateConstantParameters{ .A_ = 8069668.6184, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XYLENOOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLENO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 59619193.524, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO"), SpeciesByName(gas_phase, "BRO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 34326202.331999995, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "BENZO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD1"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOLO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.6 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD1"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD2"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD3"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARY"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BCARY"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BCARYO2VBS"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 120442815.2, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "C3H7O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.82 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 225830.27850000001, .B_ = 0.0, .C_ = -40.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "HCFC22") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 46069376.813999996, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEK"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEKO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1385092.3747999999, .B_ = 0.0, .C_ = -170.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "TOLOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TOLO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "COF2"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 12887381.2264, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOP"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPAO2"), 0.6 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPBO2"), 0.4 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 15296237.5304, .B_ = 0.0, .C_ = 410.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "E90") })
                              .SetProducts({})
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1.29e-07, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "MTERP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MTERP"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0508005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.1149005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0348005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0554005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.12780049999999998 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 379.39486788, .B_ = 0.0, .C_ = -580.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "MCO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1204428.152, .B_ = 0.0, .C_ = 500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3CL"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1180339.58896, .B_ = 0.0, .C_ = -1200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "ISOPBO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HYDRALD"), 0.87 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPNITB"), 0.08 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.92 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.92 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.05 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2649741.9343999997, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.9 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.9 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.9 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COOH"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1204428.152, .B_ = 0.0, .C_ = 500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "MACR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.12 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.24 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.65 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.88 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.33 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.14 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 903.321114, .B_ = 0.0, .C_ = -2100.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ALKOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALKO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CH3CHO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2788251.1718800003, .B_ = 0.0, .C_ = 350.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARYO2VBS"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.2202005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.20670049999999998 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0653005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.12840049999999997 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.114 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 165608.87089999998, .B_ = 0.0, .C_ = 1300.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "MTERP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NTERPO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 722656.8912, .B_ = 0.0, .C_ = 490.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N"), SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1987306.4508, .B_ = 0.0, .C_ = -3150.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HOCH2OO"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 265.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 16861994.128, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH2BR2"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 154769017.532, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CH3COCHO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 505859.82383999997, .B_ = 0.0, .C_ = 830.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "COFCL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 114420674.44, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3CCL3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 3.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 987631.08464, .B_ = 0.0, .C_ = -1520.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "HCFC141B") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COFCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 108037205.2344, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH4"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 5419926.683999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C3H7O2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C3H7OOH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOP"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOP"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPO2VBS"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 15296237.5304, .B_ = 0.0, .C_ = 410.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "TERPROD1") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERP2O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 34326202.331999995, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 873210.4102, .B_ = 0.0, .C_ = 220.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 84309970.64, .B_ = 0.0, .C_ = -470.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNO3"), SpeciesByName(gas_phase, "CH3O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NC4CHO"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NC4CH2OH"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 301107.038, .B_ = 0.0, .C_ = 400.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MDIALO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.33 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.07 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.14 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.07 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.07 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 258952.05268, .B_ = 0.0, .C_ = 1040.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "ISOPBO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.75 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYDRALD"), 0.75 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 301107.038, .B_ = 0.0, .C_ = 400.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "DMS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6624354.836, .B_ = 0.0, .C_ = -280.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HOCH2OO"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "MTERP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MTERP"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.1749305 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.5901905 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 722656.8912, .B_ = 0.0, .C_ = 490.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H6"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2H5O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4612959.82216, .B_ = 0.0, .C_ = -1020.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "GLYOXAL"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6925461.874000001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C6H5O2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C6H5OOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "MTERP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7226568.9120000005, .B_ = 0.0, .C_ = 440.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "XOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 915365.39552, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4817712.608, .B_ = 0.0, .C_ = -2060.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "IEPOX") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7828782.988, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H6"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "C2H5O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 43359413.471999995, .B_ = 0.0, .C_ = -70.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALKO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALKNIT"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 32519.560103999996, .B_ = 0.0, .C_ = 870.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BZOO"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BZALD"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 10839853.367999999, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "DMS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 114420.67444, .B_ = 0.0, .C_ = 520.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MCO3"), SpeciesByName(gas_phase, "MCO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1385092.3747999999, .B_ = 0.0, .C_ = 530.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "ENEO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2890627.5648, .B_ = 0.0, .C_ = 120.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "O2"), SpeciesByName(gas_phase, "M"), SpeciesByName(gas_phase, "O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "M"), 1.0 } })
              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                  .A_ = 217.59707599952029, .B_ = -2.4, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3854170.0864, .B_ = 0.0, .C_ = 290.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 481771.2608, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 18066422.28, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOLUENE"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TOLUENE"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TOLUO2VBS"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1023763.9292, .B_ = 0.0, .C_ = 352.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "MACR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MACRO2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MCO3"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 5781255.1296, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CFC12"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 72506574.75039999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "BRO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 11442067.444, .B_ = 0.0, .C_ = 230.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4155277.1244, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HBR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3312177.418, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARY"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 120442815.2, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MACRO2"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.22 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.47 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.53 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.22 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.53 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8430997.064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "H2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 72265689.12, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO"), SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 96354.25216, .B_ = 0.0, .C_ = -2280.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "C6H5OOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C6H5O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "HOCL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 102376.39292, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PHENO"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C6H5O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 168619.94128, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "BENZO2VBS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0097005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0034005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.1579005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0059004999999999995 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.0536005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PHENO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.7 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8430997.064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "ACBZO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C6H5O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4516605.57, .B_ = 0.0, .C_ = 290.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8430997.064, .B_ = 0.0, .C_ = 270.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "POOH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PHENO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MCO3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COOH"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COOOH"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.45 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.45 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.45 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.45 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 258952.05268, .B_ = 0.0, .C_ = 1040.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "C2H5O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "CFC11") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COFCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 124658313.73200001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "C2H5O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.7 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "C2H5OH"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 120442.8152, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPBO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 481771.2608, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "H2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1686199.4128, .B_ = 0.0, .C_ = -1800.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 963542.5216, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH4"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1475424.4862, .B_ = 0.0, .C_ = -1775.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPO2VBS"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0003005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0003005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0073005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0057005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.0623005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 350.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 28906275.648000002, .B_ = 0.0, .C_ = 250.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PHENO"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({})
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1264649.5596, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BRO"), SpeciesByName(gas_phase, "BRO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 903321.1140000001, .B_ = 0.0, .C_ = 230.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "NC4CH2OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NOA"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 42154985.32, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4456384.1624, .B_ = 0.0, .C_ = 270.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "NTERPO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPNIT"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.6 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 0.8 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2529299.1192, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "TERPROD1") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERP2O2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NTERPO2"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 602214.076, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "S"), SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1385092.3747999999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H4"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.63 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.13 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.13 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.37 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7226.568912, .B_ = 0.0, .C_ = -2630.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ST80_25") })
                              .SetProducts({})
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4.63e-07, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLONO2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 722656.8912, .B_ = 0.0, .C_ = -330.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "OCS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 12646495.595999999, .B_ = 0.0, .C_ = -2200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNO3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPNOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 481771.2608, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "MVK") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MACRO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2487144.13388, .B_ = 0.0, .C_ = 452.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPBO2"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HYDRALD"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8430997.064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H5O2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2H5OOH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H2402"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 72265689.12, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOP"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOP"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0590245 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.0250245 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1824708.65028, .B_ = 0.0, .C_ = -446.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CH3OOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.7 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3071291.7876, .B_ = 0.0, .C_ = 210.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEKOOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEKO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "MCO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3191734.6028, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2071616.42144, .B_ = 0.0, .C_ = 260.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HCFC141B") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COFCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 752767.595, .B_ = 0.0, .C_ = -1600.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BRO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1385092.3747999999, .B_ = 0.0, .C_ = 260.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ISOPNITA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.7 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.7 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.7 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HONITR"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.3 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 24088563.04, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "PO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2529299.1192, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZO2VBS"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0023005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0008005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0843005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0443005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.16210049999999998 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6022.14076, .B_ = 0.0, .C_ = -490.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CHBR3"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 3.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 278222903.112, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "MACRO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HONITR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 78287.82988, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "CF3BR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 27099633.419999998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZOOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BENZO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "TERPO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPNIT"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.32 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.04 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.8 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2529299.1192, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "S") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 39746129.016, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3COOOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 602214.076, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 19873064.508, .B_ = 0.0, .C_ = 55.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HNO3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 14453.137824, .B_ = 0.0, .C_ = 460.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HOCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 290.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "MDIALO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.83 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.17 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.35 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.17 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.17 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4516605.57, .B_ = 0.0, .C_ = 290.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H2O"), SpeciesByName(gas_phase, "F") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HF"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8430997.064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "CH3COCHO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 843099.7064, .B_ = 0.0, .C_ = -1860.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "N") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "N2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 12646495.595999999, .B_ = 0.0, .C_ = 100.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3OH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1746420.8204, .B_ = 0.0, .C_ = -345.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "XO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.1 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 301107.038, .B_ = 0.0, .C_ = 400.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "PHENOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNO3"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NC4CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8430997.064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "XYLENO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.34 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.54 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD1"), 0.06 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD2"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD3"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD4"), 0.21 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "XYLOLOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLOLO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HYAC") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1806642.2280000001, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "CH2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 361328.4456, .B_ = 0.0, .C_ = -2058.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MALO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4516605.57, .B_ = 0.0, .C_ = 290.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALKO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALKOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "TERPNIT") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 12044281.52, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HOCH2OO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2400000000000.0, .B_ = 0.0, .C_ = -7000.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BRO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 246907.77116, .B_ = 0.0, .C_ = 290.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPAO2"), SpeciesByName(gas_phase, "CH3O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MACR"), 0.31 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MVK"), 0.44 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 301107.038, .B_ = 0.0, .C_ = 400.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "H2O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 843099.7064, .B_ = 0.0, .C_ = -2000.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "EO"), SpeciesByName(gas_phase, "O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6022.14076, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N2O"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "N2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 27942733.1264, .B_ = 0.0, .C_ = 20.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HPALD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 11201181.8136, .B_ = 0.0, .C_ = 175.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HOBR"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 72265689.12, .B_ = 0.0, .C_ = -430.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH2O"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOCH2OO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 5841.4765372, .B_ = 0.0, .C_ = 625.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNO3"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NC4CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1445313.7824, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "IVOC"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "IVOCO2VBS"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8069668.6184, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2107749.266, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3CL"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 12224945.742800001, .B_ = 0.0, .C_ = -1100.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CHBR3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 3.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 541992.6684, .B_ = 0.0, .C_ = -360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "HCFC142B") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 78287829.88, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1987306.4508, .B_ = 0.0, .C_ = -115.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OCLO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 210774.9266, .B_ = 0.0, .C_ = -1370.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 21679706.735999998, .B_ = 0.0, .C_ = -375.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "C3H7O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.82 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.27 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2529299.1192, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H5OH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4155277.1244, .B_ = 0.0, .C_ = -230.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "S") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7226568.9120000005, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ISOPOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "IEPOX"), 0.6 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.6 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 9153653.9552, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "RO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1746420.8204, .B_ = 0.0, .C_ = 300.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "BZOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BZOO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 72265689.12, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MALO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.16 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.16 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.16 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 258952.05268, .B_ = 0.0, .C_ = 1040.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NH_50") })
                              .SetProducts({})
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2.31e-07, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "N2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "N2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 12947602.634, .B_ = 0.0, .C_ = 110.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "XYLENES") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLOL"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TEPOMUC"), 0.23 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BZOO"), 0.06 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLENO2"), 0.56 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.38 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 10237639.292, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 13248709.672, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "IVOCO2VBS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.1056005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.1026005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0521005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0143005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.0166005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARYO2VBS"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.1279005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.17920049999999998 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0676005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.079 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.1254005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "EO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "EO"), 0.75 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2529299.1192, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "HBR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 18066422.28, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 16861994.128, .B_ = 0.0, .C_ = 85.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH4"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4275719.9396, .B_ = 0.0, .C_ = -1270.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH2BR2"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3793948.6788000003, .B_ = 0.0, .C_ = -800.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HOCL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1806642.2280000001, .B_ = 0.0, .C_ = -500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "CH3CHO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 843099.7064, .B_ = 0.0, .C_ = -1900.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MACRO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MACROOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 481771.2608, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BR"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HBR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2890627.5648, .B_ = 0.0, .C_ = -310.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "MACRO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.47 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.53 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.53 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.22 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.22 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOP"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOP"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0033005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6323.247797999999, .B_ = 0.0, .C_ = -2000.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO"), SpeciesByName(gas_phase, "OCLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1144206.7444, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ACBZO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C6H5O2"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.4 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 258952.05268, .B_ = 0.0, .C_ = 1040.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XYLEO2VBS"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.16770049999999997 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0174005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.086 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0512005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.15980049999999998 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BZALD"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ACBZO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3553063.0484, .B_ = 0.0, .C_ = 225.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH4"), SpeciesByName(gas_phase, "F") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HF"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 96354252.16, .B_ = 0.0, .C_ = -260.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NH_5") })
                              .SetProducts({})
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2.31e-06, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPAO2"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MACR"), 0.39 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MVK"), 0.61 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8430997.064, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3COCH3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "RO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 23004577.7032, .B_ = 0.0, .C_ = -2000.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3COCH3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "RO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 80094.47210799999, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2NO2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 270996.3342, .B_ = 0.0, .C_ = 610.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7828782.988, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ALKO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MEK"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4034834.3092, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNITB"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NOA"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HONITR"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 24088563.04, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "XYLOL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLOLO2"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.63 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO"), 0.07 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 50585982.383999996, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 15657565.976, .B_ = 0.0, .C_ = 330.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MCO3"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2770184.7495999997, .B_ = 0.0, .C_ = 530.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "XYLENES") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLENES"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLEO2VBS"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 10237639.292, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HOCL"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2047527.8584, .B_ = 0.0, .C_ = -130.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CRESOL"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO2"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.73 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO"), 0.07 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 28304061.572, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CH3COOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 18969.743393999997, .B_ = 0.0, .C_ = 920.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1083985.3368, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 18066422.28, .B_ = 0.0, .C_ = -2450.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XO2"), SpeciesByName(gas_phase, "CH3CO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 782878.2988, .B_ = 0.0, .C_ = 640.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "XO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.25 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1445313.7824, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP2O2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ONITR"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.9 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.34 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.27 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.225 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.9 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD2"), 0.9 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.9 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.225 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2529299.1192, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3BR"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 855143.98792, .B_ = 0.0, .C_ = -1150.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "GLYALD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.8 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6022140.76, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NTERPO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NTERPOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "POOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PO2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "CF2CLBR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 58715872.41, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "RO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ROOH"), 0.85 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.15 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 517904.10536, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "MACRO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.47 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.22 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.53 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.22 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.53 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1445313.7824, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "NTERPO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPNIT"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.75 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1204428.152, .B_ = 0.0, .C_ = 500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MACROOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MCO3"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MACRO2"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 13850923.748000002, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MEKO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MEKOOH"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.2 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2047527.8584, .B_ = 0.0, .C_ = -1100.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZENE"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BENZENE"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BENZO2VBS"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1385092.3747999999, .B_ = 0.0, .C_ = -193.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PAN"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 24088.56304, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERPO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "MCO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3011070.38, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 13850923.748000002, .B_ = 0.0, .C_ = -200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 10237639.292, .B_ = 0.0, .C_ = 125.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "EO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "EOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3CN"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 469726.97928, .B_ = 0.0, .C_ = -1050.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1023763.9292, .B_ = 0.0, .C_ = -940.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 72265.68912, .B_ = 0.0, .C_ = -2450.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "H2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 18367529.318, .B_ = 0.0, .C_ = -2270.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "C3H6") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.12 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COOH"), 0.12 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.56 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.28 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH4"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.28 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.36 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3914.391494, .B_ = 0.0, .C_ = -1900.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C6H5O2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "ISOPBO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYDRALD"), 0.95 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.05 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1445313.7824, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N2O"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 43720741.9176, .B_ = 0.0, .C_ = 20.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "NC4CHO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NOA"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 60221407.6, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "MTERPO2VBS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0245005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.008200500000000001 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0772005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0332005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.13 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "XYLEO2VBS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0063005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0237005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0025005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.011 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.1185005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "BRO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 10237639.292, .B_ = 0.0, .C_ = 250.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "BRO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 5299483.868799999, .B_ = 0.0, .C_ = 260.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 180664.2228, .B_ = 0.0, .C_ = 460.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants(
                  { SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "M") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "M"), 1.0 } })
              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                  .A_ = 761.5897659983209, .B_ = 0.0, .C_ = 920.0, .D_ = 300.0, .E_ = 0.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants(
                  { SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "HO2"), SpeciesByName(gas_phase, "H2O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                  .A_ = 152.3179531996642, .B_ = 0.0, .C_ = 2660.0, .D_ = 300.0, .E_ = 0.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2"),
                                              SpeciesByName(gas_phase, "HO2"),
                                              SpeciesByName(gas_phase, "M"),
                                              SpeciesByName(gas_phase, "H2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "M"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 0.642096108110429, .B_ = 0.0, .C_ = 3120.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "C3H6") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NOA"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 277018.47495999996, .B_ = 0.0, .C_ = -1156.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARY"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NTERPO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 11442067.444, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARY"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 0.33 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD2"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.63 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.57 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.23 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.27 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.52 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.34 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALK"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.06 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "RO2"), 0.06 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7226.568912, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HCFC22") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 554036.9499199999, .B_ = 0.0, .C_ = -1560.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 72265689.12, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H5OOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C2H5O2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H5O2"), SpeciesByName(gas_phase, "C2H5O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.6 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "C2H5OH"), 0.4 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 40950.557168, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "N") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 30110703.8, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "XO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.25 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CFC115"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 27966821.689439997, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3CO3"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COOOH"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COOH"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.45 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.45 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 258952.05268, .B_ = 0.0, .C_ = 1040.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "PHENOL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO2"), 0.14 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO"), 0.06 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 283040.61572, .B_ = 0.0, .C_ = 1220.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCN"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 65039120.208000004, .B_ = 0.0, .C_ = 105.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNOOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NOA"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 24088563.04, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "N2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1746420.8204, .B_ = 0.0, .C_ = 220.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "HBR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 54199266.839999996, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HNO3"), SpeciesByName(gas_phase, "F") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HF"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3613284.4560000002, .B_ = 0.0, .C_ = 400.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "C3H8") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C3H7O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 5534347.35844, .B_ = 0.0, .C_ = -630.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH4"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 21077492.66, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6022140.76, .B_ = 0.0, .C_ = -3300.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCOOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 240885.6304, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BR"), SpeciesByName(gas_phase, "CH2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HBR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 10237639.292, .B_ = 0.0, .C_ = -800.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BRO"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOBR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2709963.3419999997, .B_ = 0.0, .C_ = 460.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OOH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 246907.77116, .B_ = 0.0, .C_ = 750.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNO3"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NC4CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1625978.0051999998, .B_ = 0.0, .C_ = 360.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "M") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "M"), 1.0 } })
              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                  .A_ = 100.09465495977932, .B_ = 0.0, .C_ = 720.0, .D_ = 300.0, .E_ = 0.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "MTERP") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 0.33 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD2"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.63 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.57 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.23 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.27 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.52 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.34 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALK"), 0.05 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.06 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "RO2"), 0.06 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 379.39486788, .B_ = 0.0, .C_ = -580.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "H2O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6624354.836, .B_ = 0.0, .C_ = -980.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "NOA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 403483.43091999996, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "OCS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 43359.413472, .B_ = 0.0, .C_ = -1070.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "TERPOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 19873064.508, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BIGALK"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ALKO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2107749.266, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BR"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 9635425.216, .B_ = 0.0, .C_ = -780.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZENE"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENOL"), 0.53 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BEPOMUC"), 0.12 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.65 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BENZO2"), 0.35 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1385092.3747999999, .B_ = 0.0, .C_ = -193.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "MACRO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.73 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.88 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.11 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.24 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.26 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.26 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.23 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 301107.038, .B_ = 0.0, .C_ = 400.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPAO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 481771.2608, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants(
                  { SpeciesByName(gas_phase, "H2O"), SpeciesByName(gas_phase, "H2O"), SpeciesByName(gas_phase, "SO3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2SO4"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                  .A_ = 3.082625243326537e-05, .B_ = 0.0, .C_ = 6540.0, .D_ = 300.0, .E_ = 0.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NTERPOOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NTERPO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 12044281.52, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1686199.4128, .B_ = 0.0, .C_ = 300.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPAO2"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MACR"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MVK"), 0.6 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1445313.7824, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BZOO"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BZOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "TOLUO2VBS") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.015400500000000001 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.0452005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0966005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0073005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.238 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MTERPO2VBS"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.0508005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.1149005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0348005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0554005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.12780049999999998 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 156575.65976, .B_ = 0.0, .C_ = 1300.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "MEKO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2529299.1192, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XYLENO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLENOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "DICARBO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.17 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.17 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.17 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.83 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 4516605.57, .B_ = 0.0, .C_ = 290.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP2O2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERP2OOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1987306.4508, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARY"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BCARY"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.2202005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.20670049999999998 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0653005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.12840049999999997 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.114 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 7226.568912, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "RO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 427571.99396, .B_ = 0.0, .C_ = 500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "XYLOLO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XYLOLOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BCARY"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BCARY"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.1749305 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.5901905 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 11442067.444, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O1D"), SpeciesByName(gas_phase, "CCL4") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 4.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 156997209.61319998, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "DICARBO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.07 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.07 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.07 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.33 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 258952.05268, .B_ = 0.0, .C_ = 1040.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ROOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "RO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2288413.4888, .B_ = 0.0, .C_ = 200.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "EO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 160000000000.0, .B_ = 0.0, .C_ = -4150.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOP"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MACR"), 0.3 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MVK"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.11 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.62 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.32 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.37 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.91 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.08 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "C3H6"), 0.13 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.05 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6323.247797999999, .B_ = 0.0, .C_ = -2000.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "ENEO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HONITR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 30712.917876, .B_ = 0.0, .C_ = 693.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOP"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPNO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1824708.65028, .B_ = 0.0, .C_ = -446.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH4"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 78890043.95600002, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BIGENE"), SpeciesByName(gas_phase, "NO3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.5 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 210774.9266, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "NTERPO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1445313.7824, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NH4") })
                              .SetProducts({})
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 6.34e-08, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPBO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HPALD"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1600000000.0, .B_ = 0.0, .C_ = -8300.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BENZO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BENZOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "ALKNIT") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.4 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 963542.5216, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3BR"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 8792325.5096, .B_ = 0.0, .C_ = -1040.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "TERPROD2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "RO2"), 0.15 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.68 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.8 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.65 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.2 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.7 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 20475278.584, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "XYLOLO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.17 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.51 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1565756.5976, .B_ = 0.0, .C_ = 365.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HYDRALD") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "XO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 11201181.8136, .B_ = 0.0, .C_ = 175.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP2O2"), SpeciesByName(gas_phase, "CH3O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.93 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 0.25 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.125 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.125 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.15 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1204428.152, .B_ = 0.0, .C_ = 500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CH2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3312177.418, .B_ = 0.0, .C_ = 125.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "BRONO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 11442067.444, .B_ = 0.0, .C_ = 215.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O3"), SpeciesByName(gas_phase, "MVK") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.6 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.56 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.28 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.28 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.36 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.12 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 511.8819646, .B_ = 0.0, .C_ = -1520.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL"), SpeciesByName(gas_phase, "CH2O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 48779340.156, .B_ = 0.0, .C_ = -30.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "CH3O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 301107.038, .B_ = 0.0, .C_ = -424.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOLUENE"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CRESOL"), 0.18 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TEPOMUC"), 0.1 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BZOO"), 0.07 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TOLO2"), 0.65 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.28 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1023763.9292, .B_ = 0.0, .C_ = 352.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 602214.076, .B_ = 0.0, .C_ = -1590.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "IVOCO2VBS"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 0.2381005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG1"), 0.1308005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG2"), 0.0348005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG3"), 0.0076005 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG4"), 0.0113005 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "O3") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1806642.2280000001, .B_ = 0.0, .C_ = -1500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERP2OOH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERP2O2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 13850923.748000002, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BRO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OCLO"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 572103.3722, .B_ = 0.0, .C_ = 550.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH2BR2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1204428.152, .B_ = 0.0, .C_ = -840.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "H2O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1083985.3368, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "H"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 2.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 43359413.471999995, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1083985.3368, .B_ = 0.0, .C_ = -250.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CFC113"), SpeciesByName(gas_phase, "O1D") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COFCL"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 125742299.0688, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3O2"), SpeciesByName(gas_phase, "CH3O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 11442.067444, .B_ = 0.0, .C_ = 706.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "F"), SpeciesByName(gas_phase, "H2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HF"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 84309970.64, .B_ = 0.0, .C_ = -500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NH3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 1023763.9292, .B_ = 0.0, .C_ = -710.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CHBR3"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 3.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCL"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2920738.2686, .B_ = 0.0, .C_ = -850.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLONO2"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2167970.6736, .B_ = 0.0, .C_ = -840.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPAO2"), SpeciesByName(gas_phase, "NO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPNITA"), 0.08 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.92 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MACR"), 0.36 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MVK"), 0.56 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.92 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.92 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 2649741.9343999997, .B_ = 0.0, .C_ = 180.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "H2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 9635425.216, .B_ = 0.0, .C_ = -4570.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O"), SpeciesByName(gas_phase, "HBR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3492841.6408, .B_ = 0.0, .C_ = -1500.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BIGENE"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ENEO2"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 32519560.104000002, .B_ = 0.0, .C_ = 0.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TOLO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TOLOOH"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 451660.55700000003, .B_ = 0.0, .C_ = 700.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLONO2"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::ArrheniusRateConstantParameters{
                                  .A_ = 3914391.494, .B_ = 0.0, .C_ = 135.0, .D_ = 300.0, .E_ = 0.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "TERPNIT") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jterpnit", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "O2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jo2_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3COCHO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jmgly", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "XYLENOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.34 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.54 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD1"), 0.06 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD2"), 0.2 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD3"), 0.15 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD4"), 0.21 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jxylenooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH4") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2"), 1.44 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.18 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 0.18 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.33 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 0.33 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.44 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.38 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 0.05 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch4_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "TERP2OOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.375 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.3 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.25 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.25 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jterp2ooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BRONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbrono2_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BIGALD4") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbigald4", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa3_a2") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa3_a2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BRO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbro", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "TEPOMUC") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.5 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.5 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jtepomuc", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "COF2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcof2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH2O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch2o_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BRONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbrono2_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "N"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jno", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "H2O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O1D"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jh2o_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "C3H7OOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.82 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jc3h7ooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CFC113") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COFCL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcfc113", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CCL4") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 4.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jccl4", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NOA") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jnoa", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "H2SO4") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jh2so4", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "MACR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.66 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.34 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jmacr_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NTERPOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jnterpooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "ROOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jrooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HBR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhbr", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "POOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jpooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa2_a2") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa2_a2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HYAC") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhyac", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HCFC142B") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhcfc142b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "H2O2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jh2o2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "OCLO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "joclo", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa5_a2") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa5_a2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "MACR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.34 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MCO3"), 0.66 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.34 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.34 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jmacr_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BENZOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD1"), 0.5 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbenzooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jn2o5_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "SO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "S"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jso", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "ISOPNOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "ISOPOOH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jisopnooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "ALKOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.4 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.1 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.25 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.9 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MEK"), 0.8 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jalkooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa1_a1") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa1_a1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "EOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "EO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jeooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BRCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbrcl", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "ISOPOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MVK"), 0.7 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MACR"), 0.3 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jisopooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "C6H5OOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PHENO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jc6h5ooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CFC114") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcfc114", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HO2NO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jho2no2_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "PAN") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.6 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 0.6 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.4 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 0.4 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.4 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jpan", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "O3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jo3_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CFC12") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcf2cl2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NO3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jno3_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BIGALD2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.6 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "DICARBO2"), 0.6 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbigald2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HO2NO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jho2no2_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3BR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch3br", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jclono2_b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "GLYOXAL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jglyoxal", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "O3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O1D"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jo3_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jn2o5_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HPALD") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhpald", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BZOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BZALD"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbzooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CL2O2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcl2o2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "ALKNIT") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.4 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.1 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.25 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MEK"), 0.8 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jalknit", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BIGALD") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.45 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.13 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.56 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.13 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.18 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbigald", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa5_a1") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa5_a1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "MVK") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "C3H6"), 0.7 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.7 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 0.3 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.3 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jmvk", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "SO3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jso3", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CHBR3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 3.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jchbr3", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "TERPROD2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "RO2"), 0.15 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.68 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.8 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.5 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.65 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.2 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.7 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jterprd2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "MEK") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "C2H5O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jmek", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3COCH3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jacet", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CF2CLBR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcf2clbr", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "H2O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jh2o_c", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "XYLOLOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.17 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.51 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jxylolooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "SO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jso2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "MEKOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jmekooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "PHENOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.7 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jphenooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa4_a1") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa4_a1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jno2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BIGALD3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.6 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.6 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "MDIALO2"), 0.6 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbigald3", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH2BR2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch2br2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HF") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhf", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BEPOMUC") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD1"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.5 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.5 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbepomuc", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "N2O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O1D"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "N2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jn2o", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH2O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch2o_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HCFC22") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhcfc22", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NC4CHO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jnc4cho", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HONITR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.67 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 0.33 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.33 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.33 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYALD"), 0.33 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 0.33 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.17 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.17 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhonitr", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "XOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jxooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH4") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch4_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HOBR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhobr", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NO3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jno3_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "H2O") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jh2o_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3COOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch3co3h", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "ONITR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jonitr", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "TERPOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.4 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCH3"), 0.05 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD1"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jterpooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CL2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcl2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3CHO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch3cho", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BIGALD1") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MALO2"), 0.6 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jbigald1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HNO3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhno3", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "C2H5OOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CHO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jc2h5ooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jco2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CFC115") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcfc115", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HCFC141B") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COFCL"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhcfc141b", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "TERPROD1") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "TERPROD2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jterprd1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "MPAN") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MCO3"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jmpan", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jclono2_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa3_a1") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa3_a1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa4_a2") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa4_a2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CF3BR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcf3br", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3CL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3O2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch3cl", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhcl", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "TOLOOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.6 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3COCHO"), 0.4 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD1"), 0.2 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD2"), 0.2 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "BIGALD3"), 0.2 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jtolooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3CCL3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 3.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch3ccl3", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa2_a1") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa2_a1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "GLYALD") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jglyald", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "OCS") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "S"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jocs", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "SF6") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "sink"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsf6", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HOCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jhocl", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CH3OOH") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jch3ooh", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jclo", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "COFCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "F"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcofcl", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "O2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "O"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "O1D"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jo2_a", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "H2402") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BR"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COF2"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jh2402", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "soa1_a2") })
              .SetProducts({})
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jsoa1_a2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CFC11") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 2.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "COFCL"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jcfcl3", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.5 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO"), 0.5 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 0.5 } })
              .SetRateConstant(micm::SurfaceRateConstantParameters{ .label_ = "usr_NO2_aer",
                                                                    .phase_species_ = PhaseSpeciesByName(gas_phase, "NO2"),
                                                                    .reaction_probability_ = 8e-06 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::SurfaceRateConstantParameters{ .label_ = "usr_HO2_aer",
                                                                    .phase_species_ = PhaseSpeciesByName(gas_phase, "HO2"),
                                                                    .reaction_probability_ = 0.1 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "NO3") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::SurfaceRateConstantParameters{ .label_ = "usr_NO3_aer",
                                                                    .phase_species_ = PhaseSpeciesByName(gas_phase, "NO3"),
                                                                    .reaction_probability_ = 0.002 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNITA") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_ISOPNITA_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "ISOPNITA"),
                                  .reaction_probability_ = 0.005 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ISOPNITB") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_ISOPNITB_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "ISOPNITB"),
                                  .reaction_probability_ = 0.005 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "GLYOXAL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SOAG0"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_GLYOXAL_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "GLYOXAL"),
                                  .reaction_probability_ = 0.0002 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "ONITR") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::SurfaceRateConstantParameters{ .label_ = "usr_ONITR_aer",
                                                                    .phase_species_ = PhaseSpeciesByName(gas_phase, "ONITR"),
                                                                    .reaction_probability_ = 0.005 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 2.0 } })
              .SetRateConstant(micm::SurfaceRateConstantParameters{ .label_ = "usr_N2O5_aer",
                                                                    .phase_species_ = PhaseSpeciesByName(gas_phase, "N2O5"),
                                                                    .reaction_probability_ = 0.02 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HONITR") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_HONITR_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "HONITR"),
                                  .reaction_probability_ = 0.005 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NC4CHO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_NC4CHO_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "NC4CHO"),
                                  .reaction_probability_ = 0.02 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "TERPNIT") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_TERPNIT_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "TERPNIT"),
                                  .reaction_probability_ = 0.01 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NTERPOOH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_NTERPOOH_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "NTERPOOH"),
                                  .reaction_probability_ = 0.01 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NC4CH2OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::SurfaceRateConstantParameters{
                                  .label_ = "usr_NC4CH2OH_aer",
                                  .phase_species_ = PhaseSpeciesByName(gas_phase, "NC4CH2OH"),
                                  .reaction_probability_ = 0.005 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "GLYOXAL"), 0.65 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "OH"), 0.65 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HCOOH"), 0.35 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.35 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO"), 0.35 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 1994639.8633289358,
                                                                                 .k0_B_ = 0.0,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 499837.68308,
                                                                                 .kinf_B_ = 2.0,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PAN") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CH3CO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 4.8841368205828e+17,
                                                                                 .k0_B_ = -4.1,
                                                                                 .k0_C_ = -14000.0,
                                                                                 .kinf_A_ = 1.05545e+17,
                                                                                 .kinf_B_ = -1.6,
                                                                                 .kinf_C_ = -14000.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CL2O2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLO"), 2.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 5297253446296.297,
                                                                                 .k0_B_ = -3.6,
                                                                                 .k0_C_ = -8537.0,
                                                                                 .kinf_A_ = 1712962962962963.0,
                                                                                 .kinf_B_ = -1.6,
                                                                                 .kinf_C_ = -8537.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MDIALO2"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({})
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 35178193.95325578,
                                                                                 .k0_B_ = -5.6,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5600590.9068,
                                                                                 .kinf_B_ = -1.5,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "HCN") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 2212.236939328456,
                                                                                 .k0_B_ = -1.5,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5901.6979448,
                                                                                 .kinf_B_ = 4.6,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.8,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "SO2"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 105171.9200664348,
                                                                                 .k0_B_ = -4.1,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 1023763.9292,
                                                                                 .kinf_B_ = 0.2,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 652791.2279985609,
                                                                                 .k0_B_ = -3.0,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 16861994.128,
                                                                                 .kinf_B_ = 0.0,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CLONO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 65279.12279985608,
                                                                                 .k0_B_ = -3.4,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 9033211.14,
                                                                                 .kinf_B_ = -1.9,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "C3H6") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 2901294346.66027,
                                                                                 .k0_B_ = -3.5,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 18066422.28,
                                                                                 .kinf_B_ = 0.0,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.5,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 249192041415957.1,
                                                                                 .k0_B_ = -3.0,
                                                                                 .k0_C_ = -10840.0,
                                                                                 .kinf_A_ = 275862080000000.0,
                                                                                 .kinf_B_ = 0.1,
                                                                                 .kinf_C_ = -10840.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CLO"), SpeciesByName(gas_phase, "CLO") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2O2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 6890.574073318142,
                                                                                 .k0_B_ = -3.6,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 2228192.0812,
                                                                                 .kinf_B_ = -1.6,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H4"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "EO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 31188914.226597905,
                                                                                 .k0_B_ = -3.1,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5419926.683999999,
                                                                                 .kinf_B_ = -0.85,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.48,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "O2"), SpeciesByName(gas_phase, "H") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 19221.07504662429,
                                                                                 .k0_B_ = -1.8,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 57210337.22,
                                                                                 .kinf_B_ = 0.4,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MALO2"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({})
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 35178193.95325578,
                                                                                 .k0_B_ = -5.6,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5600590.9068,
                                                                                 .kinf_B_ = -1.5,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "ACBZO2"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PBZNIT"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 35178193.95325578,
                                                                                 .k0_B_ = -5.6,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5600590.9068,
                                                                                 .kinf_B_ = -1.5,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HNO3"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 235.73016566614695,
                                                                                 .k0_B_ = 0.0,
                                                                                 .k0_C_ = 1335.0,
                                                                                 .kinf_A_ = 16.259780052,
                                                                                 .kinf_B_ = 0.0,
                                                                                 .kinf_C_ = 2199.0,
                                                                                 .Fc_ = 1.0,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H2"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 1885841.3253291757,
                                                                                 .k0_B_ = -2.4,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 132487096.72,
                                                                                 .kinf_B_ = -0.7,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MCO3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MPAN"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 35178193.95325578,
                                                                                 .k0_B_ = -5.6,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5600590.9068,
                                                                                 .kinf_B_ = -1.5,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "CH3CO3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "PAN"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 26474310.913274966,
                                                                                 .k0_B_ = -4.1,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5721033.722,
                                                                                 .kinf_B_ = -1.6,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "MPAN") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HYAC"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CH2O"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.5 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 0.5 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 2901294346.66027,
                                                                                 .k0_B_ = -3.5,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 18066422.28,
                                                                                 .kinf_B_ = 0.0,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.5,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "PBZNIT") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "ACBZO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 6.4898804328292e+17,
                                                                                 .k0_B_ = -5.6,
                                                                                 .k0_C_ = -14000.0,
                                                                                 .kinf_A_ = 1.03323e+17,
                                                                                 .kinf_B_ = -1.5,
                                                                                 .kinf_C_ = -14000.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 32639.56139992804,
                                                                                 .k0_B_ = -1.5,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 18066422.28,
                                                                                 .kinf_B_ = 0.0,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO2"), SpeciesByName(gas_phase, "HO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 68905.74073318142,
                                                                                 .k0_B_ = -3.4,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 2408856.304,
                                                                                 .kinf_B_ = -0.3,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "HO2NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 54486035447619.04,
                                                                                 .k0_B_ = -3.4,
                                                                                 .k0_C_ = -10900.0,
                                                                                 .kinf_A_ = 1904761904761904.5,
                                                                                 .kinf_B_ = -0.3,
                                                                                 .kinf_C_ = -10900.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "BRO"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRONO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 188584.13253291757,
                                                                                 .k0_B_ = -3.2,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 4155277.1244,
                                                                                 .kinf_B_ = -2.9,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "C2H4"), SpeciesByName(gas_phase, "CL") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 5802588.693320541,
                                                                                 .k0_B_ = -3.3,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 186686363.56,
                                                                                 .kinf_B_ = -1.0,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "OH") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 250236.63739944831,
                                                                                 .k0_B_ = -1.0,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 15657565.976,
                                                                                 .kinf_B_ = 0.0,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "DICARBO2"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({})
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 35178193.95325578,
                                                                                 .k0_B_ = -5.6,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 5600590.9068,
                                                                                 .kinf_B_ = -1.5,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "MPAN") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "MCO3"), 1.0 },
                                             micm::StoichSpecies{ SpeciesByName(gas_phase, "NO2"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 6.4898804328292e+17,
                                                                                 .k0_B_ = -5.6,
                                                                                 .k0_C_ = -14000.0,
                                                                                 .kinf_A_ = 1.03323e+17,
                                                                                 .kinf_B_ = -1.5,
                                                                                 .kinf_C_ = -14000.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO2"), SpeciesByName(gas_phase, "O") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "NO3"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 90665.44833313345,
                                                                                 .k0_B_ = -1.8,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 13248709.672,
                                                                                 .kinf_B_ = -0.7,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(micm::ChemicalReactionBuilder()
                              .SetReactants({ SpeciesByName(gas_phase, "NO3"), SpeciesByName(gas_phase, "NO2") })
                              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "N2O5"), 1.0 } })
                              .SetRateConstant(micm::TroeRateConstantParameters{ .k0_A_ = 870388.303998081,
                                                                                 .k0_B_ = -3.0,
                                                                                 .k0_C_ = 0.0,
                                                                                 .kinf_A_ = 963542.5216,
                                                                                 .kinf_B_ = 0.1,
                                                                                 .kinf_C_ = 0.0,
                                                                                 .Fc_ = 0.6,
                                                                                 .N_ = 1.0 })
                              .SetPhase(gas_phase)
                              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "DMS") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "SO2"), 0.5 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 0.5 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "usr_DMS_OH", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HOBR"), SpeciesByName(gas_phase, "HCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRCL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het17", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BRONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOBR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het14", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "HOCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het5", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "HOCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het10", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BRONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOBR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het3", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOCL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het2", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "BRONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOBR"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het11", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOCL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het8", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HCL"), SpeciesByName(gas_phase, "HOCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het16", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HOCL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het13", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2"), SpeciesByName(gas_phase, "HCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het9", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het1", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "HOBR"), SpeciesByName(gas_phase, "HCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "BRCL"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "H2O"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het6", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het12", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "N2O5") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 2.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het7", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2"), SpeciesByName(gas_phase, "HCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het15", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "OH"), SpeciesByName(gas_phase, "CO") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CO2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HO2"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "usr_CO_OH", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());
      processes.push_back(
          micm::ChemicalReactionBuilder()
              .SetReactants({ SpeciesByName(gas_phase, "CLONO2"), SpeciesByName(gas_phase, "HCL") })
              .SetProducts({ micm::StoichSpecies{ SpeciesByName(gas_phase, "CL2"), 1.0 },
                             micm::StoichSpecies{ SpeciesByName(gas_phase, "HNO3"), 1.0 } })
              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "het4", .scaling_factor_ = 1.0 })
              .SetPhase(gas_phase)
              .Build());

      return processes;
    }
  }  // namespace ts1

  struct Ts1
  {
    static constexpr std::string_view kName = "ts1";

    template<class Builder>
    static auto Build(Builder builder)
    {
      auto gas_phase = ts1::CreateGasPhase();
      return builder.SetSystem(micm::System(gas_phase))
          .SetReactions(ts1::CreateProcesses(gas_phase))
          .SetIgnoreUnusedSpecies(true)
          .Build();
    }

    template<class State>
    static void InitState(State& state, micm::Index num_cells)
    {
      std::vector<micm::Real> cells(num_cells);
      auto concentration = [&](const char* name, micm::Real value)
      {
        std::fill(cells.begin(), cells.end(), value);
        state[name] = cells;
      };
      auto parameter = [&](const char* label, micm::Real value)
      {
        std::fill(cells.begin(), cells.end(), value);
        state.SetCustomRateParameter(label, cells);
      };

      concentration("O2", 8.90261411);
      concentration("N2", 33.1007671);
      concentration("ALKNIT", 8.72e-10);
      concentration("BZOOH", 2.68e-11);
      concentration("C6H5OOH", 4.2e-10);
      concentration("COF2", 5.24e-11);
      concentration("COFCL", 7.92e-12);
      concentration("HF", 1.81e-11);
      concentration("F", 3.89e-28);
      concentration("BENZO2", 2.09e-13);
      concentration("BZOO", 8.81e-13);
      concentration("PAN", 4.25e-08);
      concentration("MVK", 2.61e-08);
      concentration("MACROOH", 1.21e-09);
      concentration("soa1_a1", 1e-15);
      concentration("soa1_a2", 1e-15);
      concentration("soa2_a1", 1e-15);
      concentration("soa2_a2", 1e-15);
      concentration("soa3_a1", 1e-15);
      concentration("soa3_a2", 1e-15);
      concentration("soa4_a1", 1e-15);
      concentration("soa4_a2", 1e-15);
      concentration("soa5_a1", 1e-15);
      concentration("soa5_a2", 1e-15);
      concentration("SOAG0", 3.33e-11);
      concentration("SOAG1", 2.14e-10);
      concentration("SOAG2", 1.4e-09);
      concentration("SOAG3", 4.27e-09);
      concentration("SOAG4", 2.25e-08);
      concentration("ISOPNITB", 1.05e-09);
      concentration("SO3", 1.75e-17);
      concentration("OCS", 1.96e-08);
      concentration("SO", 1.22e-20);
      concentration("S", 1.1e-26);
      concentration("H", 9.67e-18);
      concentration("MEK", 1.04e-08);
      concentration("MTERP", 1.09e-09);
      concentration("N2O5", 1e-11);
      concentration("HCN", 1.9e-08);
      concentration("SVOC", 1.98e-10);
      concentration("ISOPNO3", 5.35e-13);
      concentration("RO2", 7.13e-12);
      concentration("PHENO2", 2.54e-13);
      concentration("IVOC", 2.06e-09);
      concentration("TERPNIT", 8.28e-10);
      concentration("HNO3", 1.04e-07);
      concentration("ACBZO2", 3.54e-13);
      concentration("CH3COOOH", 7.19e-09);
      concentration("SO2", 6.15e-08);
      concentration("CH4", 7.03e-05);
      concentration("CH3CL", 2.07e-08);
      concentration("CH3CO3", 8.32e-11);
      concentration("C6H5O2", 6.23e-12);
      concentration("TERPROD1", 1.68e-09);
      concentration("HYAC", 4.88e-08);
      concentration("HPALD", 3.56e-10);
      concentration("H2O", 0.746);
      concentration("NO2", 6.37e-08);
      concentration("EOOH", 3.56e-09);
      concentration("NTERPOOH", 9.37e-11);
      concentration("CCL4", 2.96e-09);
      concentration("CF2CLBR", 1.24e-10);
      concentration("CF3BR", 1.26e-10);
      concentration("CFC11", 8.5e-09);
      concentration("CFC113", 2.66e-09);
      concentration("CFC114", 6.19e-10);
      concentration("CFC115", 3.29e-10);
      concentration("CFC12", 1.92e-08);
      concentration("CH2BR2", 4.32e-11);
      concentration("CH3BR", 2.37e-10);
      concentration("CH3CCL3", 5.55e-11);
      concentration("NO", 1.47e-08);
      concentration("BR", 3.67e-14);
      concentration("BRCL", 9.18e-19);
      concentration("BRO", 1.95e-13);
      concentration("BRONO2", 1.22e-11);
      concentration("CL", 1.07e-15);
      concentration("CL2", 1.01e-16);
      concentration("CL2O2", 1.8899999999999997e-23);
      concentration("CLO", 8.06e-14);
      concentration("CLONO2", 4.09e-11);
      concentration("HCOOH", 8.79e-09);
      concentration("HBR", 4.48e-11);
      concentration("HOBR", 1.65e-12);
      concentration("HOCL", 2.11e-12);
      concentration("BIGENE", 6.59e-10);
      concentration("C2H4", 2.51e-08);
      concentration("C2H5O2", 1.55e-12);
      concentration("CH3COCHO", 1.95e-08);
      concentration("CH3COCH3", 7.37e-08);
      concentration("O", 2.85e-14);
      concentration("OCLO", 7.68e-19);
      concentration("O1D", 1.53e-19);
      concentration("PHENO", 1e-12);
      concentration("HCFC141B", 9.29e-10);
      concentration("HCFC142B", 8.24e-10);
      concentration("HCFC22", 9.07e-09);
      concentration("DMS", 1.74e-10);
      concentration("C2H5OH", 2.75e-08);
      concentration("HCL", 8.17e-10);
      concentration("BEPOMUC", 5.4e-12);
      concentration("CHBR3", 3.87e-11);
      concentration("H2402", 1.52e-11);
      concentration("CO2", 0.0156);
      concentration("BZALD", 4.7e-10);
      concentration("BENZENE", 2.67e-09);
      concentration("C3H7O2", 2.31e-12);
      concentration("CH3O2", 3.34e-10);
      concentration("BCARY", 5.04e-12);
      concentration("BIGALD", 1.33e-11);
      concentration("BIGALD2", 6.58e-11);
      concentration("BIGALD3", 6.83e-11);
      concentration("BIGALD4", 7.34e-10);
      concentration("BIGALK", 1.9e-08);
      concentration("H2O2", 7.36e-08);
      concentration("C2H5OOH", 7.99e-11);
      concentration("C2H6", 2.93e-08);
      concentration("C3H8", 9.7e-09);
      concentration("C3H6", 2.42e-09);
      concentration("CH2O", 2.44e-07);
      concentration("CH3CN", 4.66e-09);
      concentration("C2H2", 2.13e-08);
      concentration("CH3OH", 2.74e-07);
      concentration("CH3OOH", 1.69e-08);
      concentration("CRESOL", 1.06e-10);
      concentration("ENEO2", 7.98e-12);
      concentration("MACRO2", 1.09e-10);
      concentration("ISOPAO2", 1.39e-10);
      concentration("MALO2", 1.07e-13);
      concentration("ISOPBO2", 9.08e-11);
      concentration("MCO3", 9.23e-12);
      concentration("MDIALO2", 1.61e-13);
      concentration("MEKO2", 2.81e-12);
      concentration("EO2", 4.04e-11);
      concentration("EO", 9.93e-18);
      concentration("GLYOXAL", 1.25e-08);
      concentration("MPAN", 3.11e-09);
      concentration("NC4CH2OH", 7.3e-13);
      concentration("ISOPOOH", 4.68e-09);
      concentration("GLYALD", 4.58e-08);
      concentration("HO2", 1.69e-09);
      concentration("HOCH2OO", 1.07e-13);
      concentration("H2", 1.93e-05);
      concentration("HYDRALD", 1.49e-08);
      concentration("ISOP", 1.23e-08);
      concentration("NTERPO2", 9.18e-13);
      concentration("TOLO2", 3.23e-12);
      concentration("TERP2O2", 2.15e-11);
      concentration("XYLENO2", 4.98e-12);
      concentration("TERPO2", 1.29e-11);
      concentration("XYLOLO2", 4.62e-13);
      concentration("PBZNIT", 1.78e-10);
      concentration("XYLENES", 2.82e-09);
      concentration("PO2", 1.46e-11);
      concentration("TOLUENE", 4.88e-09);
      concentration("XO2", 1.12e-10);
      concentration("XOOH", 7.93e-09);
      concentration("TERPROD2", 3.81e-09);
      concentration("MEKOOH", 6.15e-11);
      concentration("MACR", 1.06e-08);
      concentration("HONITR", 6.79e-09);
      concentration("ISOPNITA", 1.62e-09);
      concentration("ISOPNOOH", 9.16e-11);
      concentration("IEPOX", 5.84e-09);
      concentration("ONITR", 2.67e-09);
      concentration("H2SO4", 2.27e-11);
      concentration("N2O", 1.27e-05);
      concentration("NO3", 9.32e-12);
      concentration("OH", 2.52e-11);
      concentration("PHENOOH", 9.15e-12);
      concentration("PHENOL", 9.49e-11);
      concentration("XYLOL", 9.92e-11);
      concentration("ROOH", 4.4e-10);
      concentration("O3", 2.49e-06);
      concentration("TERPOOH", 2.46e-10);
      concentration("TOLOOH", 9.69e-11);
      concentration("XYLENOOH", 1.5e-10);
      concentration("XYLOLOOH", 1.31e-11);
      concentration("CO", 3.25e-06);
      concentration("CH3COOH", 8.2e-09);
      concentration("CH3CHO", 2.51e-08);
      concentration("BIGALD1", 6.47e-11);
      concentration("ALKOOH", 3.79e-10);
      concentration("DICARBO2", 1.55e-13);
      concentration("BENZOOH", 5.78e-12);
      concentration("ALKO2", 1.39e-11);
      concentration("HO2NO2", 7.99e-10);
      concentration("C3H7OOH", 7.87e-11);
      concentration("NH3", 7.59e-08);
      concentration("TERP2OOH", 5.96e-10);
      concentration("POOH", 9.24e-10);
      concentration("NOA", 5.46e-09);
      concentration("NC4CHO", 3.67e-11);
      concentration("TEPOMUC", 2.21e-10);
      concentration("NH4", 1e-15);
      concentration("SF6", 1e-15);
      concentration("E90", 1e-15);
      concentration("NH_5", 1e-15);
      concentration("NH_50", 1e-15);
      concentration("ST80_25", 1e-15);
      concentration("BENZO2VBS", 1e-15);
      concentration("ISOPO2VBS", 1e-15);
      concentration("IVOCO2VBS", 1e-15);
      concentration("MTERPO2VBS", 1e-15);
      concentration("N", 1e-15);
      concentration("TOLUO2VBS", 1e-15);
      concentration("XYLEO2VBS", 1e-15);
      concentration("BCARYO2VBS", 1e-15);
      concentration("sink", 0.0);

      parameter("jacet", 1.25e-06);
      parameter("jalknit", 5.95e-06);
      parameter("jalkooh", 5.95e-06);
      parameter("jbenzooh", 5.95e-06);
      parameter("jbepomuc", 0.00101);
      parameter("jbigald", 0.00201);
      parameter("jbigald1", 0.00141);
      parameter("jbigald2", 0.00201);
      parameter("jbigald3", 0.00201);
      parameter("jbigald4", 6.03e-05);
      parameter("jbrcl", 0.0112);
      parameter("jbro", 0.0407);
      parameter("jbrono2_a", 0.00127);
      parameter("jbrono2_b", 0.000224);
      parameter("jbzooh", 5.95e-06);
      parameter("jc2h5ooh", 5.95e-06);
      parameter("jc3h7ooh", 5.95e-06);
      parameter("jc6h5ooh", 5.95e-06);
      parameter("jccl4", 1.2799999999999999e-23);
      parameter("jcf2cl2", 2.5e-24);
      parameter("jcf2clbr", 6.95e-09);
      parameter("jcf3br", 1.45e-12);
      parameter("jcfc113", 2.6e-24);
      parameter("jcfc114", 2.09e-25);
      parameter("jcfc115", 1.2300000000000001e-26);
      parameter("jcfcl3", 1.3799999999999998e-23);
      parameter("jch2br2", 1.15e-09);
      parameter("jch2o_a", 3.47e-05);
      parameter("jch2o_b", 4.92e-05);
      parameter("jch3br", 9.35e-16);
      parameter("jch3ccl3", 1.79e-23);
      parameter("jch3cho", 7.18e-06);
      parameter("jch3cl", 4.68e-25);
      parameter("jch3co3h", 2.23e-06);
      parameter("jch3ooh", 5.95e-06);
      parameter("jch4_a", 1.0);
      parameter("jch4_b", 1.0);
      parameter("jchbr3", 1.81e-06);
      parameter("jcl2", 0.00252);
      parameter("jcl2o2", 0.00198);
      parameter("jclo", 0.00024);
      parameter("jclono2_a", 4.59e-05);
      parameter("jclono2_b", 1e-05);
      parameter("jco2", 5.28e-30);
      parameter("jcof2", 3.3e-25);
      parameter("jcofcl", 2.2999999999999997e-24);
      parameter("jeooh", 5.95e-06);
      parameter("jglyald", 6.18e-06);
      parameter("jglyoxal", 0.00014);
      parameter("jh2402", 3.68e-09);
      parameter("jh2o2", 7.98e-06);
      parameter("jh2o_a", 8.96e-31);
      parameter("jh2o_b", 1.0);
      parameter("jh2o_c", 1.0);
      parameter("jh2so4", 3.1e-10);
      parameter("jhbr", 2.5499999999999998e-23);
      parameter("jhcfc141b", 4.36e-24);
      parameter("jhcfc142b", 3.98e-26);
      parameter("jhcfc22", 8.55e-27);
      parameter("jhcl", 6.59e-25);
      parameter("jhf", 1.0);
      parameter("jhno3", 8.54e-07);
      parameter("jho2no2_a", 1.71e-06);
      parameter("jho2no2_b", 1.72e-05);
      parameter("jhobr", 0.00237);
      parameter("jhocl", 0.000304);
      parameter("jhonitr", 3.47e-05);
      parameter("jhpald", 6.03e-05);
      parameter("jhyac", 2.47e-06);
      parameter("jisopnooh", 5.95e-06);
      parameter("jisopooh", 5.95e-06);
      parameter("jmacr_a", 2.77e-06);
      parameter("jmacr_b", 2.77e-06);
      parameter("jmek", 1.25e-06);
      parameter("jmekooh", 5.95e-06);
      parameter("jmgly", 0.00014);
      parameter("jmpan", 9.69e-07);
      parameter("jmvk", 5.2e-06);
      parameter("jn2o", 9.25e-25);
      parameter("jn2o5_a", 5.54e-05);
      parameter("jn2o5_b", 3.76e-09);
      parameter("jnc4cho", 3.47e-05);
      parameter("jno", 1e-08);
      parameter("jno2", 0.0101);
      parameter("jno3_a", 0.199);
      parameter("jno3_b", 0.0124);
      parameter("jnoa", 3.47e-05);
      parameter("jnterpooh", 5.95e-06);
      parameter("jo2_a", 1.0);
      parameter("jo2_b", 1.4e-28);
      parameter("jo3_a", 5.13e-05);
      parameter("jo3_b", 0.000449);
      parameter("joclo", 0.0726);
      parameter("jocs", 6.18e-12);
      parameter("jonitr", 7.18e-06);
      parameter("jpan", 9.69e-07);
      parameter("jphenooh", 5.95e-06);
      parameter("jpooh", 5.95e-06);
      parameter("jrooh", 5.95e-06);
      parameter("jso", 2.58e-22);
      parameter("jso2", 1.18e-22);
      parameter("jso3", 2.36e-07);
      parameter("jtepomuc", 0.00101);
      parameter("jterp2ooh", 5.95e-06);
      parameter("jterpnit", 5.95e-06);
      parameter("jterpooh", 5.95e-06);
      parameter("jterprd1", 7.18e-06);
      parameter("jterprd2", 7.18e-06);
      parameter("jtolooh", 5.95e-06);
      parameter("jxooh", 5.95e-06);
      parameter("jxylenooh", 5.95e-06);
      parameter("jxylolooh", 5.95e-06);
      parameter("jsf6", 1e-20);
      parameter("jsoa1_a1", 1e-20);
      parameter("jsoa1_a2", 1e-20);
      parameter("jsoa2_a1", 1e-20);
      parameter("jsoa2_a2", 1e-20);
      parameter("jsoa3_a1", 1e-20);
      parameter("jsoa3_a2", 1e-20);
      parameter("jsoa4_a1", 1e-20);
      parameter("jsoa4_a2", 1e-20);
      parameter("jsoa5_a1", 1e-20);
      parameter("jsoa5_a2", 1e-20);
      parameter("het1", 1e-20);
      parameter("het2", 1e-20);
      parameter("het3", 1e-20);
      parameter("het4", 0.00602214076);
      parameter("het5", 0.00602214076);
      parameter("het6", 0.00602214076);
      parameter("het7", 1e-20);
      parameter("het8", 1e-20);
      parameter("het9", 0.00602214076);
      parameter("het10", 0.00602214076);
      parameter("het11", 1e-20);
      parameter("het12", 1e-20);
      parameter("het13", 1e-20);
      parameter("het14", 1e-20);
      parameter("het15", 0.00602214076);
      parameter("het16", 0.00602214076);
      parameter("het17", 0.00602214076);
      parameter("usr_NO2_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_NO2_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_HO2_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_HO2_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_NO3_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_NO3_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_ISOPNITA_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_ISOPNITA_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_ISOPNITB_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_ISOPNITB_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_GLYOXAL_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_GLYOXAL_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_ONITR_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_ONITR_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_N2O5_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_N2O5_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_HONITR_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_HONITR_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_NC4CHO_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_NC4CHO_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_TERPNIT_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_TERPNIT_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_NTERPOOH_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_NTERPOOH_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_NC4CH2OH_aer.effective radius [m]", 2.8209479177387802e-08);
      parameter("usr_NC4CH2OH_aer.particle number concentration [# m-3]", 100000000000.0);
      parameter("usr_CO_OH", 0.0);
      parameter("usr_DMS_OH", 0.0);

      for (micm::Index c = 0; c < num_cells; ++c)
      {
        state.conditions_[c].temperature_ = 287.45;
        state.conditions_[c].pressure_ = 101319.9;
        // Troe rates and the third body read air density, and nothing else
        // computes it.
        state.conditions_[c].CalculateIdealAirDensity();
      }
    }
  };
}  // namespace bench
