// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file model.cpp
 * @brief Main program for SANSMIC
 */

#ifndef VERSION_INFO
#define VERSION_INFO "local-build"
#endif

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "libsansmic.hpp"

using namespace std;

sansmic::Model::Model() {
  prefix = "temp";
  this->init_vars();
}

/**
 * @brief Create a new model for the simulation
 * @param out_prefix the prefix to use for file outputs
 */
sansmic::Model::Model(string out_prefix) {
  prefix = out_prefix;
  this->init_vars();
}

/**
 * @brief Initialize variables to zero or default values.
 */
void sansmic::Model::init_vars() {
  salt = sansmic::Salt();
  results = sansmic::Results();
  jet_model = sansmic::JetModel();

  b_use_tstfile = true;
  b_use_outfile = false;
  verbosity = 0;

  b_initialized = false;
  b_running = false;
  b_times_up = false;
  b_use_fill_table = b_use_inj_table = false;
  b_just_saved = false;

  b_first_stage = true;
  b_obi_below_roof = true;
  b_injecting = true;

  stepNum = 0;
  stageNum = 0;

  prodCell = 1;
  injCell = 1;
  obiCell = 1;
  jetPlumeCell = 1;
  noFlowCell = 1;
  startCell = 1;
  injCellBelow = 1;
  obiCellBelow = 1;
  jetPlumeCellBelow = 1;
  maxProdOrJet = 1;
  minProdOrJet = 1;

  im = 0;
  ip = 0;
  L = 0;
  m = 0;
  nEqn = 3;
  idqi = 1;
  idqf = 1;
  izbs = 1;
  i_obi = 1;
  i_obiOld;
  jplp = 1;
  jpls = 1;

  n_cells = 0;
  n_nodes = 1;
  n_wells = 1;

  z_TD0 = 0;
  h_max = -1;
  h_insol = 0.0;
  h_obi = 0.0;

  f_insol0 = 0.0;

  print_freq = 10000000000;
  dt = 0.0;
  t_tot = 0.0;
  t_wait = 0.0;
  t_last = 0;
  timet = 0.0;
  days = 0.0;
  dayOld = -1.0;

  err_cfac = 1.0;

  C_cavAve = 0;
  C_inj = 0.0;
  C_plm = 0.0;

  Q_iTot = 0.0;
  Q_fTot = 0.0;
  Q_iOld = 0.0;
  Q_fOld = 0.0;
  Q_out = 0.0;

  L_jet = 0;
  r_inj0 = 0;
  u_inj0 = 0;

  r_cavMin = 100000.0;
  r_cavMax = 0.0;

  V_stop = 0.0;
  V_plume = 0.0;
  V_insol = 0.0;
  V_ullage = 0.0;
  V_usable = 0.0;
  V_tot = 0.0;
  V_insolVent = 0.0;

  x = 0;
  u = 0;
  p1 = 0;
  p2 = 0;

  // set default model versions
  jet_model.set_version(1);
  plumeModelVer = 1;
  tempModelVer = 0;
  f_dis0 = 1.0;

  // set physical constants
  temp = 75.0;
  beta = 0.147;
  alpha = 0.09;
  D_mol = 5.03e-5;
  D_0 = 1.142e5;
  abserr = 1.0e-2;
  relerr = 1.0e-4;
}

/**
 * @brief Configure the model with the scenario provided.
 */
void sansmic::Model::configure(Scenario scenario) {
  salt = sansmic::Salt(scenario.max_brine_sg, scenario.solid_density);
  results = sansmic::Results();
  jet_model = sansmic::JetModel(scenario.jet_model_version);
  beta = scenario.diffusion_beta;
  alpha = scenario.entrainment_coeff;
  D_mol = scenario.molecular_diffusion;
  D_0 = scenario.eddy_coefficient;
  plumeModelVer = scenario.plume_model_version;
  tempModelVer = scenario.temperature_model_version;
  relerr = scenario.relative_error;
  abserr = scenario.absolute_error;
  z_TD0 = scenario.floor_depth;
  h_max = scenario.cavern_height;
  f_insol0 = scenario.fraction_insolubles;
  f_dis0 = scenario.dissolution_factor;
  n_cells = scenario.num_cells;
  n_nodes = n_cells + 1;
  n_wells = scenario.coallescing_wells;
  dx_sep = scenario.well_separation;
  h_uso = scenario.ullage_standoff;
  dataFmt = scenario.geometry_format;
  stages = scenario.stages;
  open_outfiles(false);

  if (verbosity > 0) scenario.debug_log(fileLog);

  dz = h_max / double(n_cells);
  diffCoeff = D_mol;
  C_hat = salt.get_sg_max();
  C_wall = salt.get_solid_density();
  w_hat = salt.wt_pct(C_hat, temp);

  C_tmp = vector<double>(n_nodes + 1, 0.0);
  tunc = vector<double>(n_cells + 1, -1e10);
  aa = vector<double>(n_cells + 1, 0.0);
  voldkr = vector<double>(n_nodes + 1, 0.0);
  rcscr = vector<double>(n_nodes + 1, 0.0);
  r_cav = vector<double>(n_nodes + 1, 0.0);
  r_cav0 = vector<double>(n_nodes + 1, 0.0);
  h_cav = vector<double>(n_nodes + 1, 0.0);
  C_plume = vector<double>(n_nodes + 1, 0.0);
  u_plume = vector<double>(n_nodes + 1, 0.0);
  r_plume = vector<double>(n_nodes + 1, 0.0);
  z_cav = vector<double>(n_nodes + 1, 0.0);
  p2d = vector<double>(n_nodes + 1, 0.0);
  f_dis = vector<double>(n_nodes + 1, f_dis0);
  f_disSav = vector<double>(n_nodes + 1, f_dis0);
  f_disType = vector<int>(n_nodes + 1, 0.0);
  f_dis_prt = vector<double>(n_nodes + 1, 0.0);
  f_insol = vector<double>(n_nodes + 1, f_insol0);
  tanTheta = vector<double>(n_nodes + 1, 0.0);
  cosTheta = vector<double>(n_nodes + 1, 0.0);
  depdkr = vector<double>(n_nodes + 1, 0.0);
  x_incl = vector<double>(n_nodes + 1, 0.0);
  V_injSigned = vector<double>(n_nodes + 1, 0.0);
  ss = vector<double>(n_nodes + 1, 0.0);
  phi = vector<double>(n_nodes + 1, 180.0);
  V_saltRemove = vector<double>(n_nodes + 1, 0.0);
  vec_A = vector<double>(n_nodes + 1, 0.0);
  vec_B = vector<double>(n_nodes + 1, 0.0);
  vec_C = vector<double>(n_nodes + 1, 0.0);
  vec_D = vector<double>(n_nodes + 1, 0.0);
  amd_prt = vector<double>(n_nodes + 1, 0.0);
  akd_prt = vector<double>(n_nodes + 1, 0.0);
  ca_prt = vector<double>(n_nodes + 1, 0.0);
  cb_prt = vector<double>(n_nodes + 1, 0.0);
  dC = vector<double>(n_nodes + 1, 0.0);
  dr_prt = vector<double>(n_nodes + 1, 0.0);
  C_cav = vector<double>(n_nodes + 1, C_hat);

  results.r_0 = vector<double>(n_nodes, 0.0);
  results.h_0 = vector<double>(n_nodes, 0.0);
  results.z_0 = vector<double>(n_nodes, 0.0);

  plume_rise = new PlumeRise(dz, alpha, C_cav);

  if (dataFmt == RadiusList) {
    for (int j = 1; j <= n_nodes; j++) {
      r_cav[j] = scenario.geom_radii[j];
      r_cav0[j] = scenario.geom_radii[j];
      results.r_0[j - 1] = scenario.geom_radii[j];
    }
  } else {
    throw sansmic::UNIMPLEMENTED_GEOMETRY_IDATA;
  }

  z_cav[1] = z_TD0;
  r_cav0[1] = r_cav[1];
  r_cav0[2] = r_cav[2];
  r_cav0[n_cells + 1] = r_cav[n_cells + 1];
  for (int j = 3; j <= n_cells; j++) {
    r_cav0[j] = r_cav[j];
    r_cavMin = min(r_cavMin, r_cav[j]);
    r_cavMax = max(r_cavMax, r_cav[j]);
  }
  r_cavMin = max(r_cavMin, r_cavMax / 4.0);
  V_tot = 0;
  for (int j = 1; j <= n_cells + 1; j++) {
    V_tot += _pi_ * dz * sq(r_cav[j]) * cf_to_bbl;
  }

  h_cav[1] = 0.0;
  for (int i = 1; i <= n_cells; i++) {
    im = max(i - 1, 1);
    tanTheta[i] = (r_cav[i + 1] - r_cav[im]) / (dz * double(i + 1 - im));
    theta = atan(tanTheta[i]);
    cosTheta[i] = cos(theta);
    h_cav[i + 1] = h_cav[i] + dz;
    results.h_0[i] = h_cav[i - 1] + dz;
    depdkr[i] = z_cav[1] - h_cav[i];
    results.z_0[i - 1] = z_cav[1] - h_cav[i];
  }
  depdkr[n_cells + 1] = depdkr[n_cells] - dz;
  results.z_0[n_cells] = results.z_0[n_cells - 1] - dz;
  tanTheta[n_nodes] = tanTheta[n_cells];
  V_injSigned[n_nodes] = 0.0;
  ss[n_nodes] = 0.0;
  cosTheta[n_nodes] = cosTheta[n_cells];
  phi[n_nodes] = 180.0;
  V_saltRemove[n_nodes] = 0.0;
  C_cavAve = C_hat;

  write_tst_header();  // It's nice to have this header at each stage
  save_results(results, false);
  b_initialized = true;
}

/**
 * @brief Move on to the next stage
 * @return the current stage number
 */
int sansmic::Model::init_stage(void) {
  sansmic::Stage stage = stages[stageNum];

  // save the OBI (will be 0.0 if this is the first stage)
  h_obiSave = h_obi;

  // initialize variables that have to be reset every stage;
  V_stop = -1.0e5;
  z_obi_stop = 1.0e9;
  b_injecting = true;
  Q_inj_table_name = "CONSTANT";
  b_use_inj_table = false;
  Q_fill_table_name = "CONSTANT";
  b_use_fill_table = false;

  // get the stage-defined variables
  runMode = stage.leach_mode;
  print_freq = stage.print_freq;
  b_first_stage = !(bool)stage.is_subsequent;
  t_wait = (double)stage.t_rest;
  stopCriteria = stage.stop_cond_val;
  Q_in = stage.Q_raw;
  r_tbgID = stage.r_tbgID;
  r_innerPipeOD = stage.r_tbgOD;
  r_outerPipeID = stage.r_csgID;
  r_csgOD = stage.r_csgOD;
  C_inj = stage.sg_raw;
  C_cav0 = stage.sg_init;
  h_inj = stage.h_inj;
  h_prd = stage.h_prod;
  h_obi = stage.h_obi;
  dt = stage.dt;
  t_end = stage.t_stage;
  Q_fill = stage.Q_oil;

  if (stopCriteria > 0) {
    V_stop = double(stopCriteria);
  } else if (stopCriteria < -0.0) {
    z_obi_stop = -double(stopCriteria);
  }

  // if the current depths are below zero (?) set to the
  // height of the insolubles plus the absolute height given
  // I.e., consider that these are relative heights above the
  // level of the insoluble materials not absolute heights.
  if (h_inj < 0) {
    h_inj = h_insol - h_inj;
  }
  if (h_prd < 0) {
    h_prd = h_insol - h_prd;
  }
  if (h_obi < 0) {
    h_obi = h_insol - h_obi;
  }
  if (h_uso < 0) {
    h_uso = h_insol - h_uso;
  }

  // use previous OBI if OBI is ~0
  if (h_obi < 1.0e-3) {
    h_obi = h_obiSave;
  }

  // bound concentrations
  C_inj = bound(C_inj, 1.0, C_hat);
  C_cav0 = bound(C_cav0, 1.0, C_hat);

  // get the plume radius based on the outer casing OD and inner casing ID
  r_csgOD = r_csgOD * in_to_ft;
  r_tbgID = r_tbgID * in_to_ft;
  r_inj0 = 2.0 * r_csgOD;

  // calculate the injection rates
  Q_in = Q_in * bbl_to_cf * per_d_to_per_h;
  Q_iSav = Q_in;
  Q_fill = Q_fill * bbl_to_cf * per_d_to_per_h;
  Q_fillSav = Q_fill;

  // adjust injection depth based on the jet length
  u_inj0 = jet_model.velocity(Q_in, r_tbgID);
  L_jet = jet_model.length(u_inj0);
  h_inj = max(h_inj - L_jet, 0.0);

  /**
   * b_0 = Max{l_jet, 2 r_co}
   * where r_inj0 is b_0
   */
  r_inj0 = max(L_jet, r_inj0);
  u_inj0 = jet_model.velocity(Q_in, r_inj0);  // injection point velocity

  // determine injection, production, obi cell locations
  h_inj = max(h_inj, h_insol);
  injCell = int(h_inj / dz) + 1;
  injCell = max(injCell, 1);
  prodCell = int(h_prd / dz) + 1;
  obiCell = int(h_obi / dz + 0.5) + 1;
  obiCell = min(obiCell, n_nodes);
  obiCellBelow = obiCell - 1;
  injCellBelow = max(injCell - 1, 1);

  // Begin initialization of cavern settings
  if (b_first_stage) {
    dt_min = dt;

    std::fill(C_cav.begin(), C_cav.end(), C_cav0);
    C_cav[n_nodes] = C_hat;
    C_cavAve = C_cav0;

    if (runMode == Withdrawal) {
      dt_min = (_pi_ * sq(r_cavMin) * dz) / (Q_in + 1.0e-10);
      if (dt_min < 0.01) {
        dt_min = 0.01;
      }
    }

    dt = min(dt, dt_min);
  } else {
    // initialize variables associated with the cavern sg
    for (int i = 1; i <= n_nodes; i++) {
      f_disSav[i] = f_dis[i];
    }
    // FIXME: do we really want this next line?
    if (runMode == Withdrawal && C_cav0 > 1.0) {
      std::fill(C_cav.begin(), C_cav.end(), C_cav0);
      C_cav[n_nodes] = C_hat;
      C_cavAve = C_cav0;
    }

    std::fill(V_injSigned.begin(), V_injSigned.end(), 0.0);
    std::fill(x_incl.begin(), x_incl.end(), 0.0);
    std::fill(V_saltRemove.begin(), V_saltRemove.end(), 0.0);
    std::fill(ss.begin(), ss.end(), 0.0);
    std::fill(phi.begin(), phi.end(), 180.0);
    // calculate the effective radius for coalescing caverns
    if (n_wells == 2 || n_wells == 3) {
      double p1, p2, p3, p4, p5;
      for (int i = 1; i <= n_nodes; i++) {
        p1 = sq(r_cav[i]);
        p4 = min(0.5 * dx_sep, r_cav[i]);
        p2 = sq(p4);
        p3 = ddim(p1, p2);
        p5 = ddim(r_cav[i], 0.57735 * dx_sep);
        if (n_wells == 3) {
          V_cell = 1.65399 * sq(p5) +
                   1.90986 * (p4 * sqrt(p3) + p1 * asin(p4 / r_cav[i]));
        } else {
          V_cell = p1 + 2.0 * (p4 * sqrt(p3) + p1 * asin(p4 / r_cav[i])) / _pi_;
        }
        r_cav[i] = sqrt(V_cell);
      }
      V_insol = (double)n_wells * V_insol;
    }
  }
  // End of initialization of cavern settings

  // set the prod and min timestep if withdrawal leach
  if (runMode == Withdrawal) {
    prodCell = obiCell;
    dt_min = (_pi_ * sq(r_cavMin) * dz) / (Q_in + 1.0e-10);
    // next two lines are commented out?
    //  ************************* FIXME
    // dt_min = max(dt_min, 0.01);
    // dt = min(dt, dt_min);
    // LOGGING
    // fileLog << " NOTE: THE MAXIMUM RECOMMENDED TIMESTEP IS " << dt_min
    //         << " HRS" << std::endl;
  }
  stepNum = 0;

  // temperature correction calculation
  // FIXME: the temperature is fixed at 75 degF
  for (int i = 1; i <= n_nodes; i++) {
    if (f_dis[i] > 0.92) {
      f_dis[i] =
          pow(((temp + 460.0) / 535.0), 0.75) * exp(0.01 * (temp - 75.0));
      f_disType[i] = 1;
    }
  }
  jpls = injCell + 1;

  if (n_wells > 1) {
    m_brine = brine_mass(dz);
  } else if (b_first_stage || abs(h_obiSave - h_obi) >= 0.01) {
    m_brine = brine_mass(dz);
  }

  t = 0.0;
  V_saltRemove[obiCell] = 0.0;

  // stage has been initialized
  b_running = true;  // This says the stage has been initialized
  stageNum++;        // Increment the stage number
  return stageNum;
}

/**
 * @brief End the current stage
 * @return 1 if the stage is complete
 */
int sansmic::Model::end_stage(void) {
  // Increase the time, write the log, and advance
  timet = timet + days;
  t_last = t_tot;
  h_obiSave = h_obi;
  if (verbosity > 1) {
    fileLog << "stageNum, t_last, timet, h_obi = " << stageNum << ", " << t_last
            << ", " << timet << ", " << h_obiSave << endl;
    write_log_end_stage();
  }
  return 1;
}

/**
 * @brief Run the complete model.
 */
void sansmic::Model::run_sim(void) {
  // do a complete run from initialization to finalization
  open_outfiles(false);
  for (unsigned long istage = 0; istage < stages.size(); istage++) {
    run_stage();
  }
  close_outfiles();
}

/**
 * @brief Run the next stage and all timesteps therein
 * @return 1 if the stage is complete
 */
int sansmic::Model::run_stage() {
  if (!b_running) {
    init_stage();
  }
  int status = 0;
  while (b_running) {
    status = run_step();
  }
  if (status == 0) end_stage();
  return 1;
}

/**
 * @brief Run the next time step
 * @return 1 if the stage is complete
 */
int sansmic::Model::run_step() {
  int tstep;
  if (!b_running) {
    init_stage();
  }
  tstep = leach();
  if (!b_running) {
    end_stage();
    return 1;
  }
  return 0;
}

/**
 * @brief Perform the leaching
 * @return is the stage complete
 */
int sansmic::Model::leach() {
  if (obiCell < 2) {
    std::cerr << "OBI cell = " << obiCell << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              << std::endl;
    std::cerr << "ERROR: OBI at cavern floor" << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              << std::endl;
    throw sansmic::OBI_AT_CAVERN_FLOOR;
  }

  // FIXME : injection flow table
  if (b_use_inj_table) {
    throw sansmic::UNIMPLEMENTED_FLOW_TABLES;
  }

  // FIXME: fill flow table
  if (b_use_fill_table) {
    throw sansmic::UNIMPLEMENTED_FLOW_TABLES;
  }

  // Convert back to BPD? Ugh.
  Q_inBPD = Q_in * cf_to_bbl * per_h_to_per_d;
  Q_fillBPD = Q_fill * cf_to_bbl * per_h_to_per_d;
  Q_iTot = Q_iTot + Q_in * dt * cf_to_bbl;
  Q_fTot = Q_fTot + Q_fill * dt * cf_to_bbl;

  // reinitialize temp variables for every step
  std::fill(vec_A.begin(), vec_A.end(), 0);
  std::fill(vec_B.begin(), vec_B.end(), 0);
  std::fill(vec_C.begin(), vec_C.end(), 0);
  std::fill(vec_D.begin(), vec_D.end(), 0);
  std::fill(amd_prt.begin(), amd_prt.end(), 0);
  std::fill(akd_prt.begin(), akd_prt.end(), 0);
  std::fill(ca_prt.begin(), ca_prt.end(), 0);
  std::fill(cb_prt.begin(), cb_prt.end(), 0);
  std::fill(dC.begin(), dC.end(), 0);
  std::fill(dr_prt.begin(), dr_prt.end(), 0);

  // DO PLUME MODEL EVERY 3 STEPS
  if (stepNum % 3 == 0) {
    plume(C_inj, h_inj, h_max, x, u, dz, alpha, n_nodes, r_stag, cpl);
    if (prodCell <= injCell)  // if reverse leach
    {
      jetPlumeCell = max(jetPlumeCell, jpls - 2);
    }
    jpls = jetPlumeCell;
  }

  // do a lot of bounding on variables
  jetPlumeCell = min(jetPlumeCell, obiCellBelow);
  injCell = min(injCell, obiCellBelow);
  jetPlumeCellBelow = max(jetPlumeCell - 1, 1);
  m_plume = sq(r_cav[jetPlumeCell]) * dz * C_cav[jetPlumeCell];
  vtpl = dz * double(jetPlumeCell - 1);
  dum = m_plume * (1.0 - ddim(h_insol, vtpl) / dz);
  m_plume = max(dum, 1.0);
  C_bar = C_cav[jetPlumeCell];
  V_plume = sq(r_cav[jetPlumeCell]) * dz;

  // Adjust the plume model (direct)
  if (jetPlumeCellBelow >= injCell) {
    for (int i = injCell; i <= jetPlumeCellBelow; i++) {
      r_sqd = sq(r_cav[i]) * dz;
      V_plume = V_plume + r_sqd;
      duml = (jetPlumeCell - i) * dz;
      dum = max(max(duml, r_cav[injCell] * 0.5), 0.01);
      if (prodCell > injCell) {
        if (f_dis[i] >= 0.9 * f_disSav[i]) {
          f_dis[i] = f_disSav[i] * (0.9 + 0.05 * r_cav[injCell] / dum);
          f_disType[i] = 2;
        }
      }
      // totWeightPlume is the volume of the plume * sg of the
      // plume / _pi_
      m_plume = m_plume + r_sqd * C_cav[i];
    }
    C_bar = m_plume / V_plume;
    // direct leaching velocity adjustment factor
    if (prodCell > jetPlumeCell) {
      for (int i = jetPlumeCell; i <= prodCell; i++) {
        vel = Q_in / (_pi_ * sq(r_cav[i]) + 1.0e-21);
        vel = min(vel, 4700.0);
        velfac = 1.0;
        if (vel > 100.0) {
          velfac = 0.98077 - 0.436376e-3 * vel + 0.53574e-6 * sq(vel) -
                   0.661453e-10 * vel * vel * vel;
          velfac = min(velfac, 4.0);
          f_dis[i] = f_disSav[i] * velfac;
          f_disType[i] = 3;
        }
      }
    }
  }
  V_plume = V_plume * _pi_;
  // END: ADJUST THE PLUME MODEL (direct)

  maxProdOrJet = max(prodCell, jetPlumeCell);
  minProdOrJet = min(prodCell, jetPlumeCell);
  sig = double(prodCell - jetPlumeCell) /
        dabs(double(prodCell - jetPlumeCell) + 1.0e-31);
  w_u = sig + 1.0;
  w_d = 1.0 - sig;
  for (int i = minProdOrJet; i <= maxProdOrJet; i++) {
    V_injSigned[i] = Q_in * sig;
  }

  // ADJUST THE PLUME MODEL (reverse)
  if (prodCell <= injCell) {
    double duml, dumd, dumv, p1;
    int izib;
    for (int i = injCell; i <= obiCellBelow; i++) {
      f_dis[i] = min(f_disSav[i], 1.0);
      f_disType[i] = 4;
    }
    // REVERSE LEACH ADJUSTMENT FACTOR
    if (jetPlumeCell > injCell) {
      duml = (h_obi - h_inj) * 1.15 / dz;
      izib = int(duml);
      duml = double(izib) + 0.01;
      duml = min(duml, h_obi);
      dumd = min(min(duml * dz * 0.005, 2.0), r_cav[injCell] * 0.4) *
             sqrt(Q_in) / 120.0;
      if (duml < 1.0) {
        dumd = 0.0;
      }
      izib = obiCellBelow - int(duml);
      izib = max(izib, 1);
      int ir = 0;
      int jr = 0;
      for (int i = izib; i <= obiCellBelow; i++) {
        im = max(i - 1, 1);
        ir = ir + 1;
        jr = max(izib - ir, prodCell);
        jr = min(jr, izib);
        dumv = duml - double(obiCellBelow - i);
        dum = dumv / duml;
        p1 = exp(-2.5 * dum);
        dum = min(dum, 0.999999999999999);
        dum = dumd * pow((1.0 - 0.4 * dum) * dum, 0.333);
        if (r_cav[i] > r_cav[i + 1] * 1.1 && r_cav[i] > r_cav[im]) {
          dum = 0.25 * dum;
        }
        f_dis[i] = f_disSav[i] + dum;
        f_dis[jr] = (f_dis[izib] - f_disSav[izib]) * p1 + f_disSav[jr];
        f_disType[i] = 51;
        f_disType[jr] = 52;
      }
    }
  }
  // end of adjust the plume model (reverse)

  // leach fill from production to OBI
  if (runMode >= LeachFill) {
    for (int i = prodCell; i <= obiCell; i++) {
      if (i > maxProdOrJet) {
        V_injSigned[i] = 0.0;
      }
      V_injSigned[i] = V_injSigned[i] - Q_fill;
    }
  }
  c1 = 1.0 / (dz * dz);
  c2 = 0.5 / dz;
  c3 = 2.0 * c1;
  if (prodCell <= injCell) {
    cpl = C_inj;
  }
  jplp = min(jetPlumeCell + 1, obiCell);
  cpln = C_cav[jplp];
  dt_dz2 = dt / (dz * dz);

  // calculate wall angle correction factor
  slope(jetPlumeCell);
  f_dis_prt[jetPlumeCell] = f_dis[jetPlumeCell];

  // evaluate eqn 4.1 and multiple by dis factor at jet cell
  dr_dt = salt.recession_rate(C_cav[jetPlumeCell]) * 60.0 *
          x_incl[jetPlumeCell] * f_dis[jetPlumeCell];
  dr_prt[jetPlumeCell] = dr_dt * dt;
  V_saltRemove[jetPlumeCell] =
      2.0 * dr_dt * r_cav[jetPlumeCell] * dz * dt * C_wall;
  m_saltRemove = V_saltRemove[jetPlumeCell];

  if (jetPlumeCell >= injCell) {
    for (int i = injCell; i <= jetPlumeCellBelow; i++) {
      slope(i);
      f_dis_prt[i] = f_dis[i];
      dr_dt = salt.recession_rate(C_cav[i]) * 60.0 * x_incl[i] * f_dis[i];
      dr_prt[i] = dr_dt * dt;
      V_saltRemove[i] = 2.0 * dr_dt * r_cav[i] * dz * dt * C_wall;
      m_saltRemove = m_saltRemove + V_saltRemove[i];
    }
  }

  // couple plume model
  diffCoeff = D_mol + D_0 * beta * beta;
  // if sg gradient stable or at OBI use molecular diffusion
  if (cpln <= C_cav[jetPlumeCell] || jplp == obiCell) {
    diffCoeff = D_mol;
  }

  vflo = Q_in * dt * 0.726 / (m_plume * _pi_);
  vppl = 0.726 * dt * sq(r_cav[jetPlumeCell]) * diffCoeff / (m_plume * dz);

  C_plm = (C_bar + vflo * (C_inj - C_bar) + m_saltRemove * 0.726 / m_plume +
           vppl * cpln) /
          (1.0 + vppl);

  C_plm = bound(C_plm, C_inj, C_hat);
  C_cav[jetPlumeCell] = C_plm * err_cfac;
  C_cav[jetPlumeCell] = bound(C_cav[jetPlumeCell], C_inj, C_hat);
  if (jetPlumeCell >= obiCell - 1) {
    C_cav[obiCell] = C_cav[jetPlumeCell];
  }

  for (int i = 1; i <= obiCellBelow; i++) {
    im = max(i - 1, 1);
    slope(i);
    f_dis_prt[i] = f_dis[i];

    dr_dt = salt.recession_rate(C_cav[i]) * 60.0 * f_dis[i];
    dr_prt[i] = dr_dt * dt * x_incl[i];

    volRemoved =
        _pi_ * dz * (sq(r_cav[i] + dr_dt * x_incl[i] * dt) - sq(r_cav[i]));
    if (i < minProdOrJet || i > maxProdOrJet) {
      V_injSigned[i] = 0.0;
    }
    if (i > maxProdOrJet && runMode == LeachFill) {
      V_injSigned[i] = -Q_fill;
    }
    r = max(r_cav[i], 1.0e-10);
    p5 = V_cell;
    V_cell = _pi_ * (sq(r) - sq(r_csgOD));
    V_cell = max(V_cell, 1.0e-10);
    if (i == 1) {
      p5 = V_cell;
    }
    if (i == maxProdOrJet) {
      p9 = V_cell;
    }
    p10 = V_cell;
    p1 = p5 / V_cell;
    watsol = V_cell * dz * C_cav[i];
    watsal = volRemoved * C_wall;
    w = (watsol * salt.wt_pct(C_cav[i]) + watsal) / (watsol + watsal);
    dC_dt = (salt.sg(w) - C_cav[i]) / dt;

    // diffusion
    //
    //                  (  dC  ) 1/2               2
    //  D = D_mol + D_0 ( ---- )^     (Min {r, l})^
    //                  (  dz  ) +
    //
    //                  (  dC  ) -1/4
    //  where    l ~ beta  ( ---- )^
    //                  (  dz  )
    ca = sqrt(0.5 * (ddim(C_cav[i + 1], C_cav[i]) + ddim(C_cav[i], C_cav[im])) /
              dz);
    cb = beta / (sqrt(ca) + 1.0e-10);
    p2 = diffCoeff;
    diffCoeff = D_mol + D_0 * ca * sq(min(r, cb));
    if (i == prodCell) {
      diffCoeff = diffCoeff * 1.5;
    }
    if (i == prodCell + 1) {
      diffCoeff = diffCoeff * 1.5;
    }

    if (i == 1) {
      p2 = diffCoeff;
    }

    akd_prt[i] = diffCoeff;
    amd_prt[i] = -2.0 * akd_prt[i] / r_cav[i] * 0.5 *
                 (r_cav[i + 1] - r_cav[im]) / dz * V_cell;

    ca_prt[i] = ca;
    cb_prt[i] = cb;

    // S_d = S_d from SAND2015-   ch 2.3
    //
    //
    S_d = dC_dt * V_cell * cosTheta[i] /
          (2.0 * _pi_ * r_cav[i] * diffCoeff * (C_hat - C_cav[i]) + 1.0e-7);
    S_d = max(S_d, 0.0);

    //
    ca = 2.0 * tanTheta[i] * diffCoeff / r - V_injSigned[i] / V_cell;
    p4 = ca;
    sig = ca / (dabs(ca) + 1.0e-31);
    w_u = 1.0 - sig;
    w_d = 1.0 + sig;

    // EQN 2.8 in SAND2015 ch 2.1
    cb = -(2.0 * _pi_ * r * S_d * diffCoeff / (V_cell * cosTheta[i]) +
           1.0 / dt + ss[i]);
    cc = 2.0 * _pi_ * r * diffCoeff * S_d * C_hat / (V_cell * cosTheta[i]) +
         C_cav[i] / dt + ss[i] * C_inj;

    // \bar{A_i} = D_{i-1}^n ( A_{i-1}^n / A_i^n ) / \Delta z^2
    //                 - w_u u_i^n / 2 \Delta z
    vec_A[i] = (p1 * p2 * c1) - (ca * w_u * c2);
    vec_B[i] = cb + (c2 * ca * (w_u - w_d)) - (c1 * (diffCoeff + (p2 * p1)));
    vec_C[i] = (c1 * diffCoeff) + (ca * c2 * w_d);
    vec_D[i] = -cc;

    // change the radius
    remove_salt(i, dr_dt, dt, dz);

    // account for insolubles
    r_cav[i] = sqrt(sq(r_cav[i]) + V_saltRemove[i] * f_insol[i] / (_pi_ * dz));
    if (i > maxProdOrJet) {
      p10 = p9;
    }
    m = min(i, maxProdOrJet);
    p8 = ddim(V_injSigned[m], 0.0) / p10;
    // TODO: verify this equation
    fallf = 0.5 / (1.0 + 0.00231 * p8) + 0.5 * exp(-0.002 * p8);
    fallf = min(fallf, 1.0);
    V_insol = V_insol + V_saltRemove[i] * f_insol[i] * fallf;
    V_insolVent = V_insolVent + V_saltRemove[i] * f_insol[i] * (1.0 - fallf);

    if (h_cav[i] < h_insol) {
      // THIS IS COMMENTED OUT IN THE FORTRAN CODE
      // cavRadius[i] = 0.0;
    }
  }

  if (prodCell < injCell) {
    for (int i = prodCell; i <= obiCell; i++) {
      f_dis[i] = f_disSav[i];
    }
  }
  noFlowCell = obiCell - 1;
  startCell = 1;

  // force C at plume level (internal BC)
  if (jetPlumeCell != noFlowCell) {
    vec_B[jetPlumeCell] = 1.0e20;
    vec_D[jetPlumeCell] = C_cav[jetPlumeCell] * 1.0e20;
  }

  // new plume model ; reset f_dis
  if (prodCell > jetPlumeCell) {
    for (int i = jetPlumeCell; i <= prodCell; i++) {
      f_dis[i] = f_disSav[i];
    }
  }

  // force plume cells between inj and plume to be C[jpc]
  if (jetPlumeCellBelow >= injCell) {
    for (int i = injCell; i <= jetPlumeCellBelow; i++) {
      f_dis[i] = f_disSav[i];
      vec_B[i] = 1.0e20;
      vec_D[i] = C_cav[jetPlumeCell] * 1.0e20;
    }
  }
  // end of new plume

  // No flow boundary condition at top
  if (jetPlumeCell < obiCellBelow) {
    vec_B[noFlowCell] = vec_B[noFlowCell] + vec_C[noFlowCell];
  }

  // boundary conditions for mass balance equation
  if (jetPlumeCell <= 1) {
    // set bottom cell SG to plume value
    C_tmp[1] = C_cav[jetPlumeCell];
  } else {
    V_cell = _pi_ * sq(r_cav[1]) * dz;
    V_cell = max(V_cell - V_insol, 1.0);
    if (prodCell == 1) {
      // SS looks like the \Delta C_snk without the (C_inj-C^n)*dt
      // term
      ss[1] = Q_in *
              (C_cav[1] * salt.wt_pct(C_cav[1]) -
               C_cav[2] * salt.wt_pct(C_cav[2])) *
              0.726 / (V_cell * C_cav[1] * (C_cav[1] - 1.0 + 1e-7));
    }
    ca = ddim(C_cav[2], C_cav[1]) / dz;

    // CB is like \Delta C_D in 2.5 except *ed by
    // Vol*C_1/(.726*DENSAL)
    cb = V_cell * C_cav[1] * dt * D_0 * sq(beta) * ca / (dz * 0.726 * C_wall);
    volRemoved = V_saltRemove[1] + cb;

    // CN is like C^n in 2.5 except missing \Delta C_D and divided by
    // (1+~\Delta C_snk)
    C_tmp[1] = (C_cav[1] + volRemoved * C_wall * 0.726 / (V_cell * C_cav[1]) +
                ss[1] * dt * C_inj) /
               (1.0 + ss[1] * dt);
    if (prodCell == 1) {
      C_tmp[1] = max(C_tmp[1], C_cav[2]);
    }
    C_tmp[1] = bound(C_tmp[1], C_inj, C_hat);
  }

  if (jetPlumeCell >= obiCellBelow) {
    C_tmp[obiCell] = C_cav[jetPlumeCell];
  } else {
    C_tmp[obiCell] = 0.0;
  }

  //  Solve the mass balance equation; see Appendix A in SAND2015-XXXX
  sansmic::trigad(startCell, noFlowCell, C_tmp, vec_A, vec_B, vec_C, vec_D);

  if (jetPlumeCell < obiCellBelow) {
    C_tmp[obiCell] = C_tmp[noFlowCell];
  }

  // swap the temporary concentration with the cavern concentration
  for (int i = 1; i <= obiCell; i++) {
    C_tmp[i] = bound(C_tmp[i], C_inj, C_hat);
    dum = C_cav[i];
    C_cav[i] = C_tmp[i];
    C_tmp[i] = dum;
  }

  // Percolate if unstable see SAND2015-XXXX ch 2.9
  for (int i = 2; i <= obiCellBelow; i++) {
    ip = i + 1;
    if (C_cav[i] < C_cav[ip]) {
      if (Q_in < 1.0e-3 || i <= jetPlumeCellBelow) {
        vdif = dabs(sq(r_cav[ip]) - sq(r_cav[i]));
        cold = C_cav[i];
        cnext = C_cav[ip];
        vi = sq(r_cav[i]) * C_cav[i];
        vip = sq(r_cav[ip]) * C_cav[ip];
        if (r_cav[ip] <= r_cav[i]) {
          w = (vdif * C_cav[i] * salt.wt_pct(C_cav[i]) +
               vip * salt.wt_pct(C_cav[ip])) /
              (vdif * C_cav[i] + vip);
          C_cav[i] = salt.sg(w);
          C_cav[ip] = cold;
        } else {
          w = (vdif * C_cav[ip] * salt.wt_pct(C_cav[ip]) +
               vi * salt.wt_pct(C_cav[i])) /
              (vdif * C_cav[ip] + vi);
          C_cav[i] = cnext;
          C_cav[ip] = salt.sg(w);
        }
      }
    }
    dC[i] = C_cav[i] - C_tmp[i];
  }
  dC[obiCell] = C_cav[obiCell] - C_tmp[obiCell];

  for (int i = obiCell; i <= n_nodes; i++) {
    phi[i] = 180.0;
    V_saltRemove[i] = 0.0;
    V_injSigned[i] = 0.0;
    f_dis[i] = f_disSav[i];
    if (i > obiCell) {
      C_cav[i] = C_hat;
      C_tmp[i] = 0.0;
    }
  }
  p7 = 0.0;
  m_brineNew = 0.0;
  m_brineOld = 0.0;
  m_saltRemove = 0.0;
  volRemoved = 0.0;
  V_insolRemain = V_insol;
  h_insol = 0.0;

  // *****************************************************************
  // !     CALCULATE THE NEW CAVERN MASS
  // !
  // *****************************************************************

  for (int i = 1; i <= obiCellBelow; i++) {
    V_cell = _pi_ * sq(r_cav[i]) * dz;
    volRemoved = volRemoved + V_cell;
    m_brineNew = V_cell * C_cav[i] + m_brineNew;
    p2 = V_saltRemove[i] * C_wall;
    p1 = (V_cell - V_saltRemove[i]) * C_tmp[i];
    m_brineOld = m_brineOld + p1;
    w = (p1 * salt.wt_pct(C_tmp[i]) + p2) / (p1 + p2);
    C_bar = salt.sg(w);
    p7 = p7 + V_cell - (p1 + p2) / C_bar;
    m_saltRemove = V_saltRemove[i] + m_saltRemove;
    // COMPENSATE FOR INSOLUBLES COVERING BOTTOM WALLS
    if (V_insolRemain > 0.01) {
      if (V_cell <= V_insolRemain) {
        f_dis[i] = 0.01;
        h_insol = h_insol + dz;
        V_insolRemain = V_insolRemain - V_cell;
        prodCell = max(prodCell, i + 1);
      } else {
        h_insol = h_insol + V_insolRemain * dz / V_cell;
        f_dis[i] = (V_cell - V_insolRemain) / V_cell;
        f_dis[i] = pow(f_dis[i], 0.6);
        f_dis[i] = max(f_dis[i], 0.01);
        V_insolRemain = 0.0;
      }
    }
  }

  h_inj = max(h_inj, h_insol);
  injCell = int(h_inj / dz) + 1;
  C_bar = m_brineOld / (volRemoved - m_saltRemove);
  m_saltRemove = m_saltRemove * C_wall;
  Q_out = Q_in - p7 / dt;
  w = (m_brineOld * salt.wt_pct(C_bar) + m_saltRemove) /
      (m_brineOld + m_saltRemove);
  C_bar = salt.sg(w);
  vrn = (m_brineOld + m_saltRemove) / C_bar;
  p3 = Q_in + (vrn - volRemoved) / dt;
  p5 = m_brineNew / (30.0 * Q_in * dt + m_brineNew);
  if (p5 > 0.8 || m_brineNew > 2.0e6) {
    p5 = 1.0;
  }

  Q_out = p5 * Q_out + (1.0 - p5) * p3;
  w = 0.4 + 0.6 * sq(m_brineNew) /
                ((50.0 * dt * Q_in) * (50.0 * dt * Q_in) + sq(m_brineNew));
  if (w > 0.8 || m_brineNew > 2.0e6) {
    w = 1.0;
  }
  Q_out = w * Q_out + (1.0 - w) * Q_in;
  p4 = Q_in * (C_inj * (1.0 + (C_bar - C_inj) / C_bar) / C_bar - 1.0) * 0.5;
  Q_out = Q_out + p4 * p5;

  w = salt.wt_pct(C_cav[prodCell]);
  rho_sat = w * C_cav[prodCell] / (w_hat * C_hat) * 100.0;
  if (rho_sat > 100.0) {
    rho_sat = min(rho_sat, 100.0);
  }

  izbs = obiCell;
  // move OBI if not Ordinary leaching mode
  if (runMode == LeachFill) {
    // move OBI if simultaneous LeachFill
    i_obi = int(h_obi / dz) + 1;
    dz_inc = Q_fill * dt / (_pi_ * sq(r_cav[i_obi]));
    h_obi = h_obi - min(dz_inc, dz);
    obiCell = int(h_obi / dz + 0.5) + 1;
    obiCellBelow = obiCell - 1;
  } else if (runMode != Ordinary) {
    // move OBI if Withdrawal or OilFill
    i_obiOld = obiCell;
    i_obi = int(h_obi / dz) + 1;
    dz_inc = Q_out * dt / (_pi_ * sq(r_cav[i_obi]));
    dz_inc = dabs(dz_inc);
    h_obi = h_obi + sign(1.0, Q_out) * min(dz_inc, dz);
    obiCell = int(h_obi / dz + 0.5) + 1;
    obiCell = min(obiCell, n_nodes);
    obiCellBelow = obiCell - 1;
    prodCell = obiCell;
    if (obiCell != i_obiOld) {
      C_cav[obiCell] = C_cav[obiCell - 1];
    }
    // used for delay before leaching
    if (obiCell > izbs) {
      tunc[izbs] = t + dt - 0.01 + timet * 24.0;
    }
  }

  // CALCULATE MASS BALANCE CORRECTION FACTOR
  m_brineNew = m_brineNew - V_insol * C_cav[1];
  if (runMode == Ordinary) {
    Q_fill = 0.0;
  } else if (runMode == Withdrawal) {
    Q_fill = -ddim(Q_out, 0.0);
  }
  p5 = Q_fill;
  p2 = (max(Q_out, 0.0) + Q_fill) * dt * 0.5;
  L = prodCell;
  if (prodCell > injCell) {
    L = min(prodCell + 1, obiCellBelow);
  }

  L = max(L, 1);
  p1 = m_brine + m_saltRemove + Q_in * dt * C_inj - p2 * (C_tmp[L] + C_cav[L]);
  p1 = max(p1, 0.001);
  m_brineNew = m_brineNew + (h_obi - double(izbs - 1) * dz) * _pi_ *
                                sq(r_cav[izbs]) * C_cav[obiCell];
  p7 = err_cfac;
  err_cfac = p1 / m_brineNew;
  err_cfac = 0.5 * (err_cfac + p7);
  m_brine = p1;
  Q_out = Q_out + p5;
  for (int i = 1; i <= obiCell; i++) {
    if (C_cav[i] < C_hat - 0.00001) {
      C_cav[i] = C_hat - ddim(C_hat, C_cav[i]) / sq(err_cfac);
    }
    C_cav[i] = max(C_cav[i], C_inj);
  }

  // t = t + dt;
  // Note, this is the time at the *end* of the timestep
  t = (stepNum + 1) * dt;
  t_tot = t_last + t;
  V_tot = 0.0;

  // calculate new wall angle data;
  for (int i = 1; i <= n_cells; i++) {
    im = max(i - 1, 1);
    tanTheta[i] = (r_cav[i + 1] - r_cav[im]) / (dz * double(1 + i - im));
    theta = atan(tanTheta[i]);
    V_tot = V_tot + _pi_ * sq(r_cav[i]) * dz;
    // note that no V_insol here, added to V_tot below instead
    voldkr[i + 1] = V_tot * cf_to_bbl;
    if (tunc[i] >= 0.0) {
      f_dis[i] = f_disSav[i];
    }
    cosTheta[i] = cos(theta);
  }

  V_tot = V_tot - V_insol;
  V_tot = V_tot * cf_to_bbl;
  tanTheta[n_nodes] = tanTheta[n_cells];
  cosTheta[n_nodes] = cosTheta[n_cells];

  // ULLAGE CALCULATION
  voldkr[1] = V_tot + V_insol * cf_to_bbl;
  for (int i = 2; i <= n_cells + 1; i++) {
    voldkr[i] = voldkr[1] - voldkr[i];
  }
  C_cavAve = 0.0;
  for (int i = 2; i <= max(2, obiCell - 1); i++) {
    C_cavAve = C_cavAve + C_cav[i];
  }
  C_cavAve = C_cavAve / max(1, (obiCell - 2));
  z_obi = z_cav[1] - h_obi;
  z_ullage = z_cav[1] - h_uso;  // z ullage = z TD - h ullage standoff
  V_used = xnterp_dkr(z_obi, depdkr, voldkr, n_nodes);
  V_usable = xnterp_dkr(z_ullage, depdkr, voldkr, n_nodes);
  V_ullage = V_usable - max(0.0, V_used);

  // PRINT ONE LINE SUMMARY TO *.tst ON 1 day increments
  days = t / 24.0;
  if (timet + days >= dayOld + 1.0 || (stepNum + 2) * dt > t_end) {
    write_tst_step(stageNum, b_injecting);
    Q_iOld = Q_iTot;
    Q_fOld = Q_fTot;
    dayOld = dayOld + 1;
  }

  // stop conditions
  if (stopCriteria > 0) {
    // stop when oil remaining < specified
    if (V_used <= V_stop + Q_inBPD * dt / 48.0) {
      t_end = t;
    }
  } else if (stopCriteria < 0) {
    // stop on OBI location
    if (b_injecting) {
      if (abs(z_obi - z_obi_stop) < 0.1 * dz) {
        t_end = t;
      } else if (runMode != LeachFill && z_obi < z_obi_stop) {
        t_end = t;
      } else if (z_obi > z_obi_stop) {
        t_end = t;
      }
    }
  }
  // stop on time - clear the handling variables
  b_times_up = false;
  b_obi_below_roof = true;
  b_just_saved = false;
  // save results if:
  // 1. step number is divisible by print_freq
  // 2. t == end of stage
  // 3. t will be beyond the end of the stage next time
  // 4. OBI height is above the top of the cavern
  if (((stepNum + 1) % print_freq == 0) || (h_obi > h_max)) {
    save_results(results, true);
    b_just_saved = true;
  }

  stepNum = stepNum + 1;

  if (t + dt > t_end + dt / 2.0) {
    b_times_up = true;
  } else if (h_obi <= h_max) {
    b_obi_below_roof = true;
  } else {
    h_obi = h_max - 0.1;
  }

  if (b_times_up || !b_obi_below_roof) {
    if (t_wait == 0) {
      // this happens if there is no waiting period
      // OR after the waiting period, when t_wait is
      // set to 0.
      b_running = false;
      write_tst_end_or_wo("=== END");
    } else {
      t_end = t + t_wait;
      write_tst_end_or_wo("--- WO");
      stopCriteria = 0;
      Q_in = 0.0;
      Q_fill = 0.0;
      t_wait = 0;
      b_injecting = false;
    }
  }
  if (!b_running && !b_just_saved) save_results(results, true);
  return stepNum;
}

/** @brief Calculate the brine mass
 * @param dz the height to calculate within
 * @return the mass
 */
double sansmic::Model::brine_mass(double dz) {
  // FIXME: handle insolubles above cell 0
  double m_brine = 0.0;

  for (int i = 1; i <= obiCell - 1; i++) {
    m_brine = m_brine + sq(r_cav[i]) * dz * C_cav[i];
  }
  m_brine = (m_brine + (h_obi - double(obiCell - 1) * dz) * sq(r_cav[obiCell]) *
                           C_cav[obiCell]) *
            _pi_;
  m_brine = m_brine - V_insol * C_cav[1];
  return m_brine;
}

/**
 * @brief Remove salt from the walls
 * @param i the cell to look at
 * @param recession_rate the recession rate
 * @param dt the timestep size
 * @param dz the vertical cell size
 */
void sansmic::Model::remove_salt(int i, double recession_rate, double dt,
                                 double dz) {
  double rnew;
  rnew = r_cav[i];
  if (recession_rate > 0.0) {
    rnew = r_cav[i] + recession_rate * x_incl[i] * dt;
  }
  V_saltRemove[i] = dz * _pi_ * (rnew * rnew - r_cav[i] * r_cav[i]);
  r_cav[i] = rnew;
}

/**
 * @brief Calculate the plume concentration
 */
void sansmic::Model::plume(double ci, double zi, double zb, double &x,
                           double &u, double dz, double alpha, int nmax,
                           double &r, double &cpl) {
  int neqn = 3;
  vector<double> y = vector<double>(neqn + 1, 0.0);

  double g_ft = std_gravity / foot;
  bool done = false;
  double c, dp, xout;
  int jpold;
  int iflag = 1;
  double flold = 0.0;
  double flnew;

  for (int i = 1; i <= nmax; i++) {
    C_plume[i] = 0.0;
    r_plume[i] = 0.0;
    u_plume[i] = 0.0;
  }
  C_plume[injCell - 1] = ci;
  r_plume[injCell - 1] = r_inj0;
  u_plume[injCell - 1] = u_inj0;

  // DKR bounded by 1
  c = C_cav[injCell];

  r = r_inj0;
  u = u_inj0;

  cpl = ci;
  y[1] = sq(r_inj0) * u_inj0;
  y[2] = sq(r_inj0) * sq(u_inj0);                // y[1] * u_inj0
  y[3] = sq(r_inj0) * u_inj0 * g_ft * (c - ci);  // f-star

  dp = dz * 0.1;
  x = zi;
  xout = zi + dp;

  jetPlumeCell = int(zi / dz) + 1;
  // flfac = y[1] * 3600.0;
  flold = y[1] * ci * salt.wt_pct(ci) * 3600.0;
  flnew = flold;

  // functions are integrated until (see SAND2015-XXXX Ch3):
  // u_i<1E-6; b_i>0.7r_i; or h_inj=z_{obi}
  // break if stop conditions already met
  if (u < 1.0e-6 || r >= r_cav[jetPlumeCell] * 0.7 || zi > dz * (obiCell - 1))
    done = true;

  while (!done) {
    // variables x and iflag are modified in this call
    // Equations 3.2 - 3.4 are integrated from SAND2015- Ch 3.
    plume_rise->solver.ode(y, x, xout, relerr, abserr, iflag);

    u = 0.0;
    jpold = jetPlumeCell;
    jetPlumeCell = int(x / dz) + 1;

    if (jetPlumeCell > jpold) flold = flnew;

    if (std::abs(y[1]) > 1e-10) u = y[2] / y[1];

    if (std::abs(u) > 1.0e-10) r = sqrt(dabs(y[1] / u));

    c = C_cav[jetPlumeCell];
    cpl = c - y[3] / (y[1] * g_ft + 1.0e-10);
    C_plume[jetPlumeCell] = cpl;
    r_plume[jetPlumeCell] = r;
    u_plume[jetPlumeCell] = u;

    xout = x + dp;

    // break if stop conditions met
    // 1. integration point has reached brine interface
    // 2. radius has reached 70% of cavern radius at jetPlumeCell
    // 3. jetPlumeCell within 1 cell of the brine interface
    // 4. integration has exited with flag > 3
    if (x >= zb || r >= r_cav[jetPlumeCell] * 0.7 ||
        jetPlumeCell >= obiCell - 1 || iflag > 3)
      break;
  }

  r = min(r, r_cav[jetPlumeCell]);

  return;
}

/**
 * @brief Calculate the new slope
 * @param i the cell to use
 */
void sansmic::Model::slope(int i) {
  double tanfac, xcos, tetar, teta, ratet;
  tanfac = 1.0 + tanTheta[i] * tanTheta[i];
  xcos = 1.0 / std::sqrt(tanfac);
  tetar = std::acos(xcos);

  if (tanTheta[i] >= 0.0) {
    x_incl[i] = std::sqrt(xcos);
    phi[i] = 90.0 - 57.29 * tetar;
  } else {
    phi[i] = 90.0 + 57.29 * tetar;
    teta = -57.29 * tetar;
    ratet = 1.0 + teta / 45.0;
    if (ratet <= 0.0) {
      ratet = -ratet;
      x_incl[i] = 1.0 + 0.22 * (1.0 + std::pow(ratet, 0.334));
    } else {
      x_incl[i] = 1.0 + 0.22 * (1.0 - std::pow(ratet, 0.334));
    }
  }
  return;
}

/**
 * @brief Get the single-timestep state of the model in a Results object.
 * @return Results object with single timestep of data.
 */
sansmic::Results sansmic::Model::get_current_state() {
  sansmic::Results new_results = sansmic::Results();
  new_results.r_0 = vector<double>(n_nodes, 0.0);
  new_results.h_0 = vector<double>(n_nodes, 0.0);
  new_results.z_0 = vector<double>(n_nodes, 0.0);
  for (int j = 0; j < n_nodes; j++) {
    new_results.r_0[j] = results.r_0[j];
    new_results.h_0[j] = results.h_0[j];
    new_results.z_0[j] = results.z_0[j];
  }

  save_results(new_results, false);
  return new_results;
}

/**
 * @brief Is this model currently running?
 * @return running status
 */
bool sansmic::Model::get_running_status(void) { return b_running; }

/**
 * @brief Get the current number of stages
 * @return the number of stages
 */
int sansmic::Model::get_num_stages(void) { return stages.size(); }

/**
 * @brief Get the stages object
 * @return vector of all stage definitions
 */
vector<sansmic::Stage> sansmic::Model::get_stages(void) { return stages; }

/**
 * @brief Get the compelte results object
 * @return Results
 */
sansmic::Results sansmic::Model::get_results(void) { return results; }

/**
 * @brief Get the current stage number
 * @return the stage number
 */
int sansmic::Model::get_current_stage(void) { return stageNum; }

/**
 * @brief Get the current time
 * @return the time
 */
double sansmic::Model::get_current_time(void) { return days + timet; }

/**
 * @brief Get the current cavern volume (in bbl)
 * @return the volume
 */
double sansmic::Model::get_current_volume(void) { return V_tot; }

/**
 * @brief Save results in a results object
 * @param new_results the object to add to
 * @param to_file output to files as well
 */
void sansmic::Model::save_results(sansmic::Results &new_results, bool to_file) {
  new_results.t.push_back(t_tot);
  new_results.dt.push_back(dt);
  new_results.step.push_back(stepNum);
  new_results.V_cavTot.push_back(V_tot);
  new_results.err.push_back(err_cfac);
  new_results.sg_cavAve.push_back(C_cavAve);
  new_results.V_injTot.push_back(Q_iTot);
  new_results.V_fillTot.push_back(Q_fTot);
  new_results.stage.push_back(stageNum);
  new_results.phase.push_back((int)b_injecting);
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

  // DO FILE OUTPUTS
  if (to_file) {
    write_out_timestep_summary();
  }
  p1 = 0.0;

  for (int i = 1, j = 1; j <= n_nodes; j++) {
    i = n_nodes + 1 - j;
    p1 = p1 + _pi_ * dz * sq(r_cav[i]) * cf_to_bbl;
    p2 = V_injSigned[i] * cf_to_bbl * per_h_to_per_d;
    if (to_file) {
      write_out_timestep_per_cell(i, p1, p2);
    }
  }
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

  if (to_file) {
    write_out_timestep_totals(Q_outBPD);
  }
  new_results.Q_out.push_back(Q_outBPD);
  new_results.sg_out.push_back(C_cav[prodCell]);

  p1 = V_insol * cf_to_bbl;
  p2 = V_insolVent * cf_to_bbl;
  if (to_file) {
    write_out_timestep_insolubles(p1, p2);
  }
  new_results.V_insolTot.push_back(p1);
  new_results.V_insolVent.push_back(p2);
  new_results.h_insol.push_back(h_insol);
  new_results.z_insol.push_back(z_cav[1] - h_insol);
  new_results.z_obi.push_back(z_cav[1] - h_obi);

  p1 = volRemoved + _pi_ * sq(r_cav[izbs]) * (h_obi - double(izbs - 1) * dz) -
       V_insol;
  p1 = p1 * cf_to_bbl;
  if (to_file) {
    write_out_timestep_removed(p1);
  }
}

/**
 * @brief Choose whether the .TST file should be written
 * @param use_file the choice
 */
void sansmic::Model::generate_tst_file(bool use_file) {
  b_use_tstfile = use_file;
}

/**
 * @brief Choose whether the .OUT file should be written
 * @param use_file the choice
 */
void sansmic::Model::generate_out_file(bool use_file) {
  b_use_outfile = use_file;
}

/**
 * @brief Set the verbosity output level for cout/cerr
 * @param verb the verbosity level
 */
void sansmic::Model::set_verbosity_level(int verb) { verbosity = verb; }

/**
 * @brief Open output files and initialize the headers.
 * @param append whether the files should be opened with append, by default
 * false
 * @return success code
 */
int sansmic::Model::open_outfiles(bool append) {
  if (!this->fileOut.is_open() && b_use_outfile) {
    if (append) {
      this->fileOut.open(prefix + ".out", std::ios::app);
    } else {
      this->fileOut.open(prefix + ".out", std::ios::out);
      write_header(this->fileOut);
    }
  }
  if (!this->fileOut.is_open() && b_use_outfile) {
    std::cerr << "Error writing file " << prefix << ".out - exiting" << endl;
    return sansmic::INVALID_OUT_FILE;
  }
  if (!this->fileLog.is_open()) {
    if (!append) {
      this->fileLog.open(prefix + ".log", std::ios::out);
      write_header(this->fileLog);
    } else {
      this->fileLog.open(prefix + ".log", std::ios::app);
    }
  }
  if (!this->fileLog.is_open()) {
    std::cerr << "Error writing file " << prefix << ".log - exiting" << endl;
    return sansmic::INVALID_LOG_FILE;
  }
  if (!this->fileTst.is_open() && b_use_tstfile) {
    if (append) {
      this->fileTst.open(prefix + ".tst", std::ios::app);
    } else {
      this->fileTst.open(prefix + ".tst", std::ios::out);
      write_header(this->fileTst);
    }
  }
  if (!this->fileTst.is_open() && b_use_tstfile) {
    std::cerr << "Error writing file " << prefix << ".tst - exiting" << endl;
    return sansmic::INVALID_TST_FILE;
  }
  return 0;
}

/**
 * @brief Close the file handles for output files
 */
void sansmic::Model::close_outfiles(void) {
  if (this->fileOut.is_open() && b_use_outfile) this->fileOut.close();
  if (this->fileLog.is_open()) this->fileLog.close();
  if (this->fileTst.is_open() && b_use_tstfile) this->fileTst.close();
}

/**
 * @brief Add header information to a file
 * @param sout the file to write to
 */
void sansmic::Model::write_header(ofstream &sout) {
  sout << "# Generated by sansmic v" << VERSION_INFO
       << " under the BSD-3-Clause license." << endl;
  sout << "# [sansmic (c) 2024 NTESS, "
          "https://github.com/SandiaLabs/sansmic/blob/main/LICENSE]"
       << endl;
}

/**
 * @brief Initialize the TST file
 */
void sansmic::Model::write_tst_header(void) {
  if (!b_use_tstfile) return;
  fileTst << "File= " << prefix << endl;
  fileTst << "#" << std::setw(12) << "t"  // time in days
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
  fileTst << " #" << std::setw(11) << "(d)" << std::setw(13) << "(bbl)"
          << std::setw(13) << "(:1)" << std::setw(13) << "(:1.kg/L)"
          << std::setw(13) << "(:1.kg/L)" << std::setw(13) << "(bbl)"
          << std::setw(13) << "(ft)" << std::setw(13) << "(ft)" << std::setw(13)
          << "(bbl)" << std::setw(12) << "(bbl)" << std::setw(13) << "(bbl)"
          << std::setw(13) << "(bbl/d)" << std::setw(13) << "(bbl)"
          << std::setw(13) << "(bbl/d)" << std::setw(13) << "(bbl)" << endl;
}

/**
 * @brief Write data to the TST file
 * @param stage the stage number
 * @param inject whether this is injection or workover
 */
void sansmic::Model::write_tst_step(int stage, bool inject) {
  if (!b_use_tstfile) return;
  fileTst << std::scientific;
  fileTst << std::setprecision(4);
  fileTst << std::setw(13) << (days + timet) << std::setw(13) << V_tot
          << std::setw(13) << err_cfac << std::setw(13) << C_cav[prodCell]
          << std::setw(13) << C_cavAve << std::setw(13) << V_insol * cf_to_bbl
          << std::setw(13) << z_cav[1] - h_insol << std::setw(13)
          << z_cav[1] - h_obi << std::setw(13) << V_insolVent * cf_to_bbl
          << std::setw(12) << V_ullage << std::setw(13) << V_usable
          << std::setw(13) << Q_iTot - Q_iOld << std::setw(13) << Q_iTot
          << std::setw(13) << Q_fTot - Q_fOld << std::setw(13) << Q_fTot
          << endl;
}

/**
 * @brief Write end of phase to the TST file
 * @param text what to use as a prefix
 */
void sansmic::Model::write_tst_end_or_wo(const char *text) {
  if (!b_use_tstfile) return;
  fileTst
      << "                                                     "
         "                                                                   "
         "                                                 "
      << text << std::setw(2) << stageNum << endl;
}

/**
 * @brief Write data to the LOG file
 */
void sansmic::Model::write_log_end_stage(void) {
  fileLog << "END STAGE ----------------------------------" << endl;
}

/**
 * @brief Write data to the OUT file
 */
void sansmic::Model::write_out_timestep_summary(void) {
  if (!b_use_outfile) return;
  fileOut << endl;
  fileOut << std::setprecision(4);
  fileOut << std::scientific;
  fileOut << " TIME=" << setw(11) << days << " DAYS (" << setw(11)
          << days * 24.0 << " HRS)"
          << "     "
          << "DT=" << setw(11) << dt << " HOURS"
          << "   START TIME=" << setw(11) << timet << " DAYS (" << setw(11)
          << timet * 24.0 << " HRS)" << endl;

  if (b_injecting) {
    fileOut << " INJECTION PHASE " << setw(5) << stageNum << endl;
  } else {
    fileOut << " WORKOVER PHASE " << setw(5) << stageNum << endl;
  }
  fileOut << endl;

  fileOut << "  I(INJ), I(PRD), I(OBI), I(PLM)= " << std::setw(12) << injCell
          << std::setw(12) << prodCell << std::setw(12) << obiCell
          << std::setw(12) << jetPlumeCell << endl;
  fileOut << "  H(INJ), H(PRD), H(OBI), H(PLM)= " << std::setw(12)
          << h_cav[injCell] << std::setw(12) << h_cav[prodCell] << std::setw(12)
          << h_cav[obiCell] << std::setw(12) << h_cav[jetPlumeCell] << endl;
  fileOut << "  Ljet(ft), RO(ft), UO(ft/s)    = " << std::setw(12) << L_jet
          << std::setw(12) << r_inj0 << std::setw(12) << u_inj0 << endl;
  fileOut << endl;
  fileOut << "     "        // cell index number
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
  fileOut << endl;
}

/**
 * @brief Write data to the OUT file
 * @param i cell index
 * @param p1 flow rate
 * @param p2 volume
 */
void sansmic::Model::write_out_timestep_per_cell(int i, double p1, double p2) {
  if (!b_use_outfile) return;
  fileOut << " " << setw(4) << i << setw(11) << h_cav[i] << setw(11) << r_cav[i]
          << setw(11) << r_cav[i] - r_cav0[i] << setw(11) << C_cav[i]
          << setw(11) << phi[i] << setw(11) << p2 << setw(11) << p1 << setw(11)
          << f_dis_prt[i] << setw(5) << f_disType[i] << setw(11) << x_incl[i]
          << setw(11) << amd_prt[i] * cf_to_bbl * per_h_to_per_d << setw(11)
          << akd_prt[i] << setw(11) << ca_prt[i] << setw(11) << C_tmp[i]
          << setw(11) << C_cav[i] << setw(11) << dC[i] << setw(11) << dr_prt[i]
          << setw(11) << C_plume[i] << setw(11) << u_plume[i] << setw(11)
          << r_plume[i] << endl;
}

/**
 * @brief Write data to the OUT file
 * @param qo output flow rate
 */
void sansmic::Model::write_out_timestep_totals(double qo) {
  if (!b_use_outfile) return;
  fileOut << endl;
  fileOut << " TIME=" << setw(11) << days << " DAYS (" << setw(11)
          << days * 24.0 << " HRS)"
          << "     "
          << "DT=" << setw(11) << dt << " HOURS"
          << "   START TIME=" << setw(11) << timet << " DAYS (" << setw(11)
          << timet * 24.0 << " HRS)" << endl;

  if (b_injecting) {
    fileOut << " INJECTION PHASE " << setw(5) << stageNum << endl;
  } else {
    fileOut << " WORKOVER PHASE " << setw(5) << stageNum << endl;
  }
  fileOut << endl;
  fileOut << " TOTAL VOLUME           =" << setw(11) << V_tot << " BBLS "
          << endl;
  fileOut << " BRINE OUT              =" << setw(11) << qo << " BBLS/DAY"
          << endl;
  fileOut << " OUTLET SPECIFIC GRAVITY=" << setw(11) << C_cav[prodCell] << endl;
}

/**
 * @brief Write data to the OUT file
 * @param p1 insolubles deposited
 * @param p2 insolubles vented
 */
void sansmic::Model::write_out_timestep_insolubles(double p1, double p2) {
  if (!b_use_outfile) return;
  fileOut << " VOLUME OF INSOLUBLES   =" << setw(11) << p1 << " BBLS " << endl;
  fileOut << " INSOL LEVEL            =" << setw(11) << h_insol << " FT"
          << endl;
  fileOut << " BLANKET LEVEL          =" << setw(11) << h_obi << " FT" << endl;
  fileOut << " VOL OF INS VENTED      =" << setw(11) << p2 << " BBLS" << endl;
}

/**
 * @brief Write data to the OUT file
 * @param p1 total brine volume
 */
void sansmic::Model::write_out_timestep_removed(double p1) {
  if (!b_use_outfile) return;
  fileOut << " BRINE VOLUME           =" << setw(11) << p1 << " BBLS " << endl;
  fileOut << " ULLAGE                 =" << setw(11) << V_ullage << " BBLS"
          << endl;
  fileOut << " USEABLE VOLUME         =" << setw(11) << V_usable << " BBLS"
          << endl;
  fileOut << " CFAC                   =" << setw(11) << err_cfac << endl;
  fileOut << "----------------------------------------------------------"
          << "----------------------------------------------------------"
          << endl;
}
