// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file scenario.cpp
 * @brief Scenario options definitions
 */

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "libsansmic.hpp"

/**
 * @brief Blank options with default version and coefficient values.
 */
sansmic::Scenario::Scenario() {
  title = "";
  num_cells = -1;
  cavern_height = -1.0;
  floor_depth = 10000.0;
  geometry_format = RadiusList;
  geom_radii.clear();
  geom_depths.clear();
  geom_volumes.clear();
  fraction_insolubles = -1.0;
  ullage_standoff = 0.0;
  jet_model_version = 1;
  plume_model_version = 1;
  temperature_model_version = 0;
  coallescing_wells = 1;
  well_separation = 0.0;
  entrainment_coeff = 0.09;
  diffusion_beta = 0.147;
  molecular_diffusion = 5.03e-5;
  eddy_coefficient = 1.142e5;
  relative_error = 1.0e-4;
  absolute_error = 1.0e-2;
  dissolution_factor = 1.0;
  max_brine_sg = 1.2019;
  solid_density = 2.16;
}

/**
 * @brief Add a stage to the simulation.
 * @param stage the stage definition struct
 * @return the current number of stages
 */
int sansmic::Scenario::add_stage(Stage stage) {
  stages.push_back(stage);
  return stages.size();
}

/**
 * @brief Add a debug output to a file
 * @param fout the open steam to print to
 */
void sansmic::Scenario::debug_log(ofstream &fout) {
  fout << "- level: DEBUG" << endl;
  fout << "  file: \"scenario.cpp:--\"" << endl;
  fout << "  funcName: \"sansmic::Scenario::debug_log\"" << endl;
  fout << "  message: Data received for scenario"
       << std::endl
       << "  scenario:" << std::endl;
  fout << "    title: '" << title << "'" << std::endl;
  fout << "    comments: '''" << comments << "'''" << std::endl;
  fout << "    num-cells: " << num_cells << std::endl;
  fout << "    cavern-height: " << cavern_height << std::endl;
  fout << "    floor-depth: " << floor_depth << std::endl;
  fout << "    ullage-standoff: " << ullage_standoff << std::endl;
  fout << "    geometry-format: " << geometry_format << std::endl;
  fout << "    fraction-insolubles: " << fraction_insolubles << std::endl;
  fout << "    geom-radii: [";
  for (int i = 0; i < geom_radii.size(); i++) fout << geom_radii[i] << ", ";
  fout << "]" << std::endl;
  fout << "    geom-depths: [";
  for (int i = 0; i < geom_depths.size(); i++) fout << geom_depths[i] << ", ";
  fout << "]" << std::endl;
  fout << "    geom-volumes: [";
  for (int i = 0; i < geom_volumes.size(); i++) fout << geom_volumes[i] << ", ";
  fout << "]" << std::endl;
  fout << "    advanced:" << std::endl;
  fout << "      coallescing-wells = " << coallescing_wells << std::endl;
  fout << "      well-separation = " << well_separation << std::endl;
  fout << "      jet-model-version = " << jet_model_version << std::endl;
  fout << "      plume-model-version = " << plume_model_version << std::endl;
  fout << "      temperature-model-version = " << temperature_model_version
       << std::endl;
  fout << "      entrainment-coeff = " << entrainment_coeff << std::endl;
  fout << "      diffusion-beta-coeff = " << diffusion_beta << std::endl;
  fout << "      molecular-diffusion-coeff = " << molecular_diffusion
       << std::endl;
  fout << "      eddy-coeff = " << eddy_coefficient << std::endl;
  fout << "      ode-relerr = " << relative_error << std::endl;
  fout << "      ode-abserr = " << absolute_error << std::endl;
  fout << "      dissolution-factor = " << dissolution_factor << std::endl;
  fout << "      salt-max-brine-sg = " << max_brine_sg << std::endl;
  fout << "      salt-solid-density = " << solid_density << std::endl;
  fout << "    stages:" << std::endl;
  for (int i = 0; i < stages.size(); i++) stages[i].debug_log(fout);
}
