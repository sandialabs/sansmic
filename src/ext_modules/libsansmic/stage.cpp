// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file stage.cpp
 * @brief Stage definition
 */

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "libsansmic.hpp"

/**
 * @brief Stage constructor with default values.
 */
sansmic::Stage::Stage() {
  title = "";
  leach_mode = Withdrawal;
  print_freq = -1;
  is_subsequent = 0;
  t_rest = -1;
  stop_cond_val = 0;
  h_inj = -1;
  h_prod = -1;
  h_obi = -1;
  Q_raw = 0;
  r_tbgID = -1;
  r_tbgOD = -1;
  r_csgID = -1;
  r_csgOD = -1;
  Q_oil = 0;
  sg_raw = 1.0;
  sg_init = 1.2019;
  dt = 0.1;
}

/**
 * @brief Output debug data to a stream.
 * @param fout Open stream to output to
 */
void sansmic::Stage::debug_log(ofstream &fout, int stageNum) {
  fout << "- level: DEBUG" << endl;
  fout << "  file: \"stage.cpp:--\"" << endl;
  fout << "  funcName: \"sansmic::Stage::debug_log\"" << endl;
  fout << "  message: Data received for stage number " << stageNum << std::endl;
  fout << "  stage:" << endl;
  fout << "    title = '" << title << "'" << endl;
  fout << "    leach-mode = " << leach_mode << endl;
  fout << "    print-freq = " << print_freq << endl;
  fout << "    is-subsequent = " << is_subsequent << endl;
  fout << "    t-rest = " << t_rest << endl;
  fout << "    stop-cond-val = " << stop_cond_val << endl;
  fout << "    h-inj = " << h_inj << endl;
  fout << "    h-prod = " << h_prod << endl;
  fout << "    h-obi = " << h_obi << endl;
  fout << "    Q-raw = " << Q_raw << endl;
  fout << "    r-tbgID = " << r_tbgID << endl;
  fout << "    r-tbgOD = " << r_tbgOD << endl;
  fout << "    r-csgID = " << r_csgID << endl;
  fout << "    r-csgOD = " << r_csgOD << endl;
  fout << "    sg-raw = " << sg_raw << endl;
  fout << "    sg-init = " << sg_init << endl;
  fout << "    dt = " << dt << endl;
  fout << "    t-stage = " << t_stage << endl;
  fout << "    Q-oil = " << Q_oil << endl;
}
