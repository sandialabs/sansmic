// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file basemodel.cpp
 * @brief Base sansmic model class.
 */

#ifndef VERSION_INFO
#define VERSION_INFO "local-build"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "libsansmic.hpp"

/**
 * @brief Is this model currently running?
 * @return running status
 */
bool sansmic::BaseModel::get_is_running(void) { return b_is_running; }

/**
 * @brief Get the current number of stages
 * @return the number of stages
 */
int sansmic::BaseModel::get_num_stages(void) { return stages.size(); }

/**
 * @brief Get the stages object
 * @return vector of all stage definitions
 */
vector<sansmic::Stage> sansmic::BaseModel::get_stages(void) { return stages; }

/**
 * @brief Get the current stage number
 * @return the stage number
 */
int sansmic::BaseModel::get_current_stage(void) { return stageNum; }

/**
 * @brief Get the current time
 * @return the time
 */
double sansmic::BaseModel::get_current_time(void) { return days + timet; }

/**
 * @brief Get the current cavern volume (in bbl)
 * @return the volume
 */
double sansmic::BaseModel::get_current_volume(void) { return V_tot; }

/**
 * @brief Choose whether the .TST file should be written
 * @param use_file the choice
 */
void sansmic::BaseModel::set_use_tstfile(bool use_file) {
  b_use_tstfile = use_file;
}

/**
 * @brief Choose whether the .OUT file should be written
 * @param use_file the choice
 */
void sansmic::BaseModel::set_use_outfile(bool use_file) {
  b_use_outfile = use_file;
}

/**
 * @brief Set the verbosity output level for cout/cerr
 * @param verb the verbosity level
 */
void sansmic::BaseModel::set_verbosity_level(int verb) { verbosity = verb; }

/**
 * @brief Get the compelte results object
 * @return Results
 */
sansmic::Results sansmic::BaseModel::get_results(void) { return results; }

/**
 * @brief Get the single-timestep state of the model in a Results object.
 * @return Results object with single timestep of data.
 */
sansmic::Results sansmic::BaseModel::get_current_state() {
  sansmic::Results new_results = sansmic::Results();
  new_results.r_0 = vector<double>(n_nodes, 0.0);
  new_results.h_0 = vector<double>(n_nodes, 0.0);
  new_results.z_0 = vector<double>(n_nodes, 0.0);
  for (int j = 0; j < n_nodes; j++) {
    new_results.r_0[j] = results.r_0[j];
    new_results.h_0[j] = results.h_0[j];
    new_results.z_0[j] = results.z_0[j];
  }

  save_results(new_results);
  return new_results;
}

/**
 * @brief Save results in a results object
 * @param new_results the object to add to
 * @param to_file output to files as well
 */
void sansmic::BaseModel::save_results(sansmic::Results &new_results) {
  double p1, p2;

  new_results.t.push_back(t_tot);
  new_results.dt.push_back(dt);
  new_results.step.push_back(stepNum);
  new_results.V_cavTot.push_back(V_tot);
  new_results.err.push_back(err_cfac);
  new_results.sg_cavAve.push_back(C_cavAve);
  new_results.V_injTot.push_back(Q_iTot);
  new_results.V_fillTot.push_back(Q_fTot);
  new_results.stage.push_back(stageNum);
  new_results.phase.push_back((int)b_is_injecting);
  new_results.injCell.push_back(injCell);
  new_results.prodCell.push_back(prodCell);
  new_results.obiCell.push_back(obiCell);
  new_results.plmCell.push_back(jetPlumeCell);
  new_results.z_plm.push_back(z_cav[1] - h_cav[jetPlumeCell]);
  new_results.z_inj.push_back(z_cav[1] - h_cav[injCell]);
  new_results.z_prod.push_back(z_cav[1] - h_cav[prodCell]);
  new_results.l_jet.push_back(L_jet);
  new_results.r_jet.push_back(r_inj0);
  new_results.u_jet.push_back(u_inj0);

  vector<double> _r_cav, _dr_cav, _sg, _theta, _Q_inj, _V, _f_dis, _f_flag,
      _xincl, _amd, _D_coeff, _dC_dz, _C_old, _C_new, _dC_dt, _dr_dt, _C_plm,
      _u_plm, _r_plm;
  for (int i = 1; i <= n_nodes; i++) {
    p1 = _pi_ * dz * sq(r_cav[i]) * cf_to_bbl;
    p2 = V_injSigned[i] * cf_to_bbl * per_h_to_per_d;
    _r_cav.push_back(r_cav[i]);
    _dr_cav.push_back(r_cav[i] - r_cav0[i]);
    _sg.push_back(C_cav[i]);
    _theta.push_back(phi[i]);
    _Q_inj.push_back(p2);
    _V.push_back(p1);
    _f_dis.push_back(f_dis_prt[i]);
    _f_flag.push_back(f_disType[i]);
    _xincl.push_back(x_incl[i]);
    _amd.push_back(amd_prt[i] * cf_to_bbl * per_h_to_per_d);
    _D_coeff.push_back(akd_prt[i]);
    _dC_dz.push_back(ca_prt[i]);
    _C_old.push_back(C_tmp[i]);
    _C_new.push_back(C_cav[i]);
    _dC_dt.push_back(dC[i]);
    _dr_dt.push_back(dr_prt[i]);
    _C_plm.push_back(C_plume[i]);
    _u_plm.push_back(u_plume[i]);
    _r_plm.push_back(r_plume[i]);
  }
  new_results.r_cav.push_back(_r_cav);
  new_results.dr_cav.push_back(_dr_cav);
  new_results.sg.push_back(_sg);
  new_results.theta.push_back(_theta);
  new_results.Q_inj.push_back(_Q_inj);
  new_results.V.push_back(_V);
  new_results.f_dis.push_back(_f_dis);
  new_results.f_flag.push_back(_f_flag);
  new_results.xincl.push_back(_xincl);
  new_results.amd.push_back(_amd);
  new_results.D_coeff.push_back(_D_coeff);
  new_results.dC_dz.push_back(_dC_dz);
  new_results.C_old.push_back(_C_old);
  new_results.C_new.push_back(_C_new);
  new_results.dC_dt.push_back(_dC_dt);
  new_results.dr_dt.push_back(_dr_dt);
  new_results.C_plm.push_back(_C_plm);
  new_results.u_plm.push_back(_u_plm);
  new_results.r_plm.push_back(_r_plm);

  slctim = days + timet;
  for (int i = 1; i <= n_nodes; i++) {
    if (h_cav[i] < h_insol) {
      rcscr[i] = 0.0;
    } else {
      rcscr[i] = r_cav[i];
    }
  }
  Q_outBPD = Q_out * cf_to_bbl * per_h_to_per_d;  // 4.27387;

  new_results.Q_out.push_back(Q_outBPD);
  new_results.sg_out.push_back(C_cav[prodCell]);

  p1 = V_insol * cf_to_bbl;
  p2 = V_insolVent * cf_to_bbl;
  new_results.V_insolTot.push_back(p1);
  new_results.V_insolVent.push_back(p2);
  new_results.h_insol.push_back(h_insol);
  new_results.z_insol.push_back(z_cav[1] - h_insol);
  new_results.z_obi.push_back(z_cav[1] - h_obi);
}

/**
 * @brief Add header information to a file
 * @param sout the file to write to
 */
void sansmic::BaseModel::write_version_info(ofstream &sout) {
  sout << "# Generated by sansmic v" << VERSION_INFO
       << " under the BSD-3-Clause license." << endl;
  sout << "# [sansmic (c) 2024 NTESS, "
          "https://github.com/SandiaLabs/sansmic/blob/main/LICENSE]"
       << endl;
}

/**
 * @brief Output a stage initialization summary to a log file
 * @param sout the output stream to write to
 */
void sansmic::BaseModel::write_simulation_init(ofstream &sout) {
  if (verbosity > 0) {
    sout << std::fixed;
    if (verbosity > 1) {
      sout << "- message: Setup simulation" << std::endl;
      sout << "  data:" << std::endl;
    } else {
      sout << "- Setup simulation:" << std::endl;
    }
    sout << std::setprecision(0);
    sout << "    initial-cavern-volume_bbl:  " << std::setw(8) << V_tot
         << std::endl;
    sout << std::setprecision(2);
    sout << "    initial-floor-depth_ftMD:   " << std::setw(11) << z_cav[1]
         << std::endl;
    sout << "    model-extent_ft:            " << std::setw(11) << h_max
         << std::endl;
    sout << std::setprecision(0);
    sout << "    number-of-nodes:            " << std::setw(8) << n_nodes
         << std::endl;
    sout << std::setprecision(2);
    sout << "    cell-height_ft:             " << std::setw(11) << dz
         << std::endl;
  }
}

/**
 * @brief Output a stage initialization summary to a log file
 * @param sout the output stream to write to
 */
void sansmic::BaseModel::write_stage_init(ofstream &sout) {
  if (verbosity > 1) {
    sout << std::fixed;
    sout << "- message: Setup stage " << stageNum << std::endl;
    sout << "  data:" << std::endl;
    sout << std::setprecision(2);
    sout << "    interface-depth_ftMD:       " << std::setw(11)
         << z_cav[1] - h_obi << std::endl;
    sout << "    production-depth_ftMD:      " << std::setw(11)
         << z_cav[1] - h_prd << std::endl;
    sout << "    injection-depth_ftMD:       " << std::setw(11)
         << z_cav[1] - h_inj << std::endl;
    sout << std::setprecision(4);
    sout << "    injected-water-conc_sg:     " << std::setw(13) << C_inj
         << std::endl;
  }
}

/**
 * @brief Output a stage completion summary to a log file
 * @param sout the output stream to write to
 */
void sansmic::BaseModel::write_stage_summary(ofstream &sout) {
  if (verbosity > 1) {
    sout << std::fixed;
    sout << "- message: End of stage" << std::endl;
    sout << "  data:" << std::endl;
    sout << std::setprecision(2);
    sout << "    cumul-sim-time_h:           " << std::setw(11) << t_tot
         << std::endl;
    sout << "    interface-level_ftMD:       " << std::setw(11)
         << z_cav[1] - h_obi << std::endl;
    sout << "    new-floor-depth_ftMD:       " << std::setw(11)
         << z_cav[1] - h_insol << std::endl;
    sout << std::setprecision(4);
    sout << "    average-brine-conc_sg:      " << std::setw(13) << C_cavAve
         << std::endl;
    sout << std::setprecision(0);
    sout << "    cumul-H2O-injected_bbl:     " << std::setw(8) << Q_iTot
         << std::endl;
    sout << "    cumul-product-in-out_bbl:  " << std::setw(9) << Q_fTot
         << std::endl;
    sout << "    cavern-volume_bbl:          " << std::setw(8) << V_tot
         << std::endl;
  }
}

/**
 * @brief Output a stage completion summary to a log file
 * @param sout the output stream to write to
 */
void sansmic::BaseModel::write_simulation_summary(ofstream &sout) {
  double r_max_final = 0;
  double dr_max_final = 0;
  if (verbosity > 0) {
    r_max_final =
        *std::max_element(results.r_cav[results.r_cav.size() - 1].begin(),
                          results.r_cav[results.r_cav.size() - 1].end());
    dr_max_final =
        *std::max_element(results.dr_cav[results.dr_cav.size() - 1].begin(),
                          results.dr_cav[results.dr_cav.size() - 1].end());
    sout << std::fixed;
    if (verbosity > 1) {
      sout << "- message: End of simulation" << std::endl;
      sout << "  data:" << std::endl;
    } else {
      sout << "- End of simulation:" << std::endl;
    }
    sout << std::setprecision(0);
    sout << "    number-of-stages-completed: " << std::setw(8) << stageNum
         << std::endl;
    sout << std::setprecision(2);
    sout << "    total-simulated-time_h:     " << std::setw(11) << t_tot
         << std::endl;
    sout << std::setprecision(0);
    sout << "    final-cavern-volume_bbl:    " << std::setw(8) << V_tot
         << std::endl;
    sout << std::setprecision(2);
    sout << "    final-interface-level_ftMD: " << std::setw(11)
         << z_cav[1] - h_obi << std::endl;
    sout << "    final-floor-depth_ft:       " << std::setw(11)
         << z_cav[1] - h_insol << std::endl;
    sout << std::setprecision(4);
    sout << "    final-brine-conc_sg:        " << std::setw(13) << C_cavAve
         << std::endl;
    sout << std::setprecision(0);
    sout << "    total-H2O-injected_bbl:     " << std::setw(8) << Q_iTot
         << std::endl;
    sout << "    total-product-in-out_bbl:  " << std::setw(9) << Q_fTot
         << std::endl;
    sout << std::setprecision(2);
    sout << "    final-max-radius_ft:        " << std::setw(11) << r_max_final
         << endl;
    sout << "    max-change-in-radius_ft:    " << std::setw(11) << dr_max_final
         << endl;
  }
}

/**
 * @brief Initialize the TST file
 */
void sansmic::BaseModel::write_daily_column_header(ofstream &sout) {
  if (b_use_tstfile) {
    sout << "File= " << prefix << endl;
    sout << "#" << std::setw(12) << "t"  // time in days
         << std::setw(13) << "V_tot"     // total volume
         << std::setw(13)
         << "err_conv"  // measure of stability of the mass balance solver
         << std::setw(13) << "sg_out"  // specific gravity of produced brine
         << std::setw(13) << "sg_ave"  // average specific gravity of the
                                       // brine within the cavern
         << std::setw(13)
         << "V_insol"  // total volume, in bbl, of insolubes create that has
                       // accumulated at the bottom
         << std::setw(13) << "D_insol"  // depth of the top of insoluble level
         << std::setw(13) << "D_OBI"    // depth of OBI
         << std::setw(13)
         << "V_insolVnt"  // volume of insolubles brought to surface entrained
                          // in the brine, in bbl
         << std::setw(12) << "V_ullage"  // volume of cavern available for oil
                                         // above ullage ref depth
         << std::setw(13)
         << "V_usable"  // volume of cavern above the ullage ref point
         << std::setw(13) << "Q_inj"   // current injection of brine in bbl/d
         << std::setw(13) << "V_inj"   // total volume injected
         << std::setw(13) << "Q_fill"  // rate of oil injection
         << std::setw(13) << "V_fill"  // total volume oil injected
         << endl;
    sout << " #" << std::setw(11) << "(d)" << std::setw(13) << "(bbl)"
         << std::setw(13) << "(:1)" << std::setw(13) << "(:1.kg/L)"
         << std::setw(13) << "(:1.kg/L)" << std::setw(13) << "(bbl)"
         << std::setw(13) << "(ft)" << std::setw(13) << "(ft)" << std::setw(13)
         << "(bbl)" << std::setw(12) << "(bbl)" << std::setw(13) << "(bbl)"
         << std::setw(13) << "(bbl/d)" << std::setw(13) << "(bbl)"
         << std::setw(13) << "(bbl/d)" << std::setw(13) << "(bbl)" << endl;
  }
}

/**
 * @brief Write data to the TST file
 * @param stage the stage number
 * @param inject whether this is injection or workover
 */
void sansmic::BaseModel::write_daily_summary(ofstream &sout, int stage,
                                             bool inject) {
  if (b_use_tstfile) {
    sout << std::scientific;
    sout << std::setprecision(4);
    sout << std::setw(13) << (days + timet) << std::setw(13) << V_tot
         << std::setw(13) << err_cfac << std::setw(13) << C_cav[prodCell]
         << std::setw(13) << C_cavAve << std::setw(13) << V_insol * cf_to_bbl
         << std::setw(13) << z_cav[1] - h_insol << std::setw(13)
         << z_cav[1] - h_obi << std::setw(13) << V_insolVent * cf_to_bbl
         << std::setw(12) << V_ullage << std::setw(13) << V_usable
         << std::setw(13) << Q_iTot - Q_iOld << std::setw(13) << Q_iTot
         << std::setw(13) << Q_fTot - Q_fOld << std::setw(13) << Q_fTot << endl;
  }
}

/**
 * @brief Write end of phase to the TST file
 * @param text what to use as a prefix
 */
void sansmic::BaseModel::write_daily_end_of_phase(ofstream &sout,
                                                  const char *text) {
  if (b_use_tstfile) {
    sout
        << "                                                     "
           "                                                                   "
           "                                                 "
        << text << std::setw(2) << stageNum << endl;
  }
}

void sansmic::BaseModel::write_detailed_results(ofstream &sout) {
  double p1 = 0.0;
  double p2 = 0.0;
  if (b_use_outfile) {
    sout << endl;
    sout << std::setprecision(4);
    sout << std::scientific;
    sout << " TIME=" << setw(11) << days << " DAYS (" << setw(11) << days * 24.0
         << " HRS)"
         << "     "
         << "DT=" << setw(11) << dt << " HOURS"
         << "   START TIME=" << setw(11) << timet << " DAYS (" << setw(11)
         << timet * 24.0 << " HRS)" << endl;

    if (b_is_injecting) {
      sout << " INJECTION PHASE " << setw(5) << stageNum << endl;
    } else {
      sout << " WORKOVER PHASE " << setw(5) << stageNum << endl;
    }
    sout << endl;

    sout << "  I(INJ), I(PRD), I(OBI), I(PLM)= " << std::setw(12) << injCell
         << std::setw(12) << prodCell << std::setw(12) << obiCell
         << std::setw(12) << jetPlumeCell << endl;
    sout << "  H(INJ), H(PRD), H(OBI), H(PLM)= " << std::setw(12)
         << h_cav[injCell] << std::setw(12) << h_cav[prodCell] << std::setw(12)
         << h_cav[obiCell] << std::setw(12) << h_cav[jetPlumeCell] << endl;
    sout << "  Ljet(ft), RO(ft), UO(ft/s)    = " << std::setw(12) << L_jet
         << std::setw(12) << r_inj0 << std::setw(12) << u_inj0 << endl;
    sout << endl;
    sout << "     "        // cell index number
         << "   HEIGHT  "  // height above initial cavern floor
         << "  RADIUS   "  // radius, in feet
         << "   d(RAD)  "  // cumulative change in radius, in feet
         << " BRINE S.G."  // sg of brine; oil set to 1.2019
         << " WALL ANGLE"  // wall angle, in degrees
         << " FLOW RATE "  // brine flow rate in bbl/d
         << "   VOLUME  "  // cumulative volume of cavern, from roof
         << "  DISFAC   "  // dissolution adjustment factor
         << "     "        // adjustemnt factor type flag
         << "  XINCL    "  // wall inclination dissolution factor
         << "   AMD     "  // internal debugging variable
         << "     AKD   "  // diffusion coefficient
         << "     dc/dz "  // change in specific gravity from cell to cell
         << "     Cold  "  // previous calvulated value of specific grav
         << "     Cnew  "  // current calculated value of sg in cell
         << "     dC    "  // Cnew - Cold
         << "    dR     "  // wall recession
         << "  Cplm     "  // brine concentration of plume
         << "  Uplm     "  // velocity of plume
         << " Rplm      "  // radius of plume
         << endl;
    sout << endl;
    for (int i = 1, j = 1; j <= n_nodes; j++) {
      i = n_nodes + 1 - j;
      p1 = p1 + _pi_ * dz * sq(r_cav[i]) * cf_to_bbl;
      p2 = V_injSigned[i] * cf_to_bbl * per_h_to_per_d;
      sout << " " << setw(4) << i << setw(11) << h_cav[i] << setw(11)
           << r_cav[i] << setw(11) << r_cav[i] - r_cav0[i] << setw(11)
           << C_cav[i] << setw(11) << phi[i] << setw(11) << p2 << setw(11) << p1
           << setw(11) << f_dis_prt[i] << setw(5) << f_disType[i] << setw(11)
           << x_incl[i] << setw(11) << amd_prt[i] * cf_to_bbl * per_h_to_per_d
           << setw(11) << akd_prt[i] << setw(11) << ca_prt[i] << setw(11)
           << C_tmp[i] << setw(11) << C_cav[i] << setw(11) << dC[i] << setw(11)
           << dr_prt[i] << setw(11) << C_plume[i] << setw(11) << u_plume[i]
           << setw(11) << r_plume[i] << endl;
    }
    Q_outBPD = Q_out * cf_to_bbl * per_h_to_per_d;  // 4.27387;
    sout << endl;
    sout << " TIME=" << setw(11) << days << " DAYS (" << setw(11) << days * 24.0
         << " HRS)"
         << "     "
         << "DT=" << setw(11) << dt << " HOURS"
         << "   START TIME=" << setw(11) << timet << " DAYS (" << setw(11)
         << timet * 24.0 << " HRS)" << endl;

    if (b_is_injecting) {
      sout << " INJECTION PHASE " << setw(5) << stageNum << endl;
    } else {
      sout << " WORKOVER PHASE " << setw(5) << stageNum << endl;
    }
    sout << endl;
    sout << " TOTAL VOLUME           =" << setw(11) << V_tot << " BBLS "
         << endl;
    sout << " BRINE OUT              =" << setw(11) << Q_outBPD << " BBLS/DAY"
         << endl;
    sout << " OUTLET SPECIFIC GRAVITY=" << setw(11) << C_cav[prodCell] << endl;
    p1 = V_insol * cf_to_bbl;
    p2 = V_insolVent * cf_to_bbl;
    sout << " VOLUME OF INSOLUBLES   =" << setw(11) << p1 << " BBLS " << endl;
    sout << " INSOL LEVEL            =" << setw(11) << h_insol << " FT" << endl;
    sout << " BLANKET LEVEL          =" << setw(11) << h_obi << " FT" << endl;
    sout << " VOL OF INS VENTED      =" << setw(11) << p2 << " BBLS" << endl;
    p1 = volRemoved + _pi_ * sq(r_cav[izbs]) * (h_obi - double(izbs - 1) * dz) -
         V_insol;
    p1 = p1 * cf_to_bbl;
    sout << " BRINE VOLUME           =" << setw(11) << p1 << " BBLS " << endl;
    sout << " ULLAGE                 =" << setw(11) << V_ullage << " BBLS"
         << endl;
    sout << " USEABLE VOLUME         =" << setw(11) << V_usable << " BBLS"
         << endl;
    sout << " CFAC                   =" << setw(11) << err_cfac << endl;
    sout << "----------------------------------------------------------"
         << "----------------------------------------------------------"
         << endl;
  }
}
