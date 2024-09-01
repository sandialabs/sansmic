// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file smoke.cpp
 * @brief Smoke implementations for the sansmic.hpp header file.
 */

#include "sansmic.hpp"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

std::string sansmic::error::get_error_text(int err) {
  return "error text";
}

int
sansmic::Model::open_outfiles (bool append)
{
  if (!this->fileOut.is_open ())
    {
      if (append)
        {
          this->fileOut.open (prefix + ".out", std::ios::app);
        }
      else
        {
          this->fileOut.open (prefix + ".out", std::ios::out);
          write_header (this->fileOut);
        }
    }
  if (!this->fileOut.is_open ())
    {
      std::cerr << "Error writing file " << prefix << ".out - exiting" << endl;
      return sansmic::error::INVALID_OUT_FILE;
    }
  if (!this->fileLog.is_open ())
    {

      if (!append)
        {
          this->fileLog.open (prefix + ".log", std::ios::out);
          write_header (this->fileLog);
        }
      else
        {
          this->fileLog.open (prefix + ".log", std::ios::app);
        }
    }
  if (!this->fileLog.is_open ())
    {
      std::cerr << "Error writing file " << prefix << ".log - exiting" << endl;
      return sansmic::error::INVALID_LOG_FILE;
    }
  if (!this->fileTst.is_open ())
    {
      if (append)
        {
          this->fileTst.open (prefix + ".tst", std::ios::app);
        }
      else
        {
          this->fileTst.open (prefix + ".tst", std::ios::out);
          write_header (this->fileTst);
        }
    }
  if (!this->fileTst.is_open ())
    {
      std::cerr << "Error writing file " << prefix << ".tst - exiting" << endl;
      return sansmic::error::INVALID_TST_FILE;
    }
  return 0;
}

void
sansmic::Model::close_outfiles (void)
{
  if (this->fileOut.is_open ())
    this->fileOut.close ();
  if (this->fileLog.is_open ())
    this->fileLog.close ();
  if (this->fileTst.is_open ())
    this->fileTst.close ();
}

void
sansmic::Model::write_header (ofstream &sout)
{
  sout << "[[write_header]]" << std::endl;
  sout << "header = 'smoke test'" << std::endl;
}

void
sansmic::Model::write_tst_header (void)
{
  fileTst << "[[write_tst_header]]" << std::endl;
  fileTst << "header = 'tst header'" << std::endl;
}

void
sansmic::Model::write_tst_step (int stage, bool inject)
{
  fileTst << "[[write_tst_step]]" << std::endl;
  fileTst << "stage = " << stage << std::endl;
  fileTst << "inject = " << inject << std::endl;
}

void
sansmic::Model::write_tst_end_or_wo (const char *text)
{
  fileTst << "[[write_tst_end_or_wo]]" << std::endl;
  fileTst << "text = " << text << std::endl;
}

void
sansmic::Model::write_log_end_stage (void)
{
  fileLog << "[[write_log_end_stage]]" << std::endl;
  fileLog << "end-stage = true" << std::endl;
}

void
sansmic::Model::write_out_timestep_summary (void)
{
  fileOut << "[[write_out_timestep_summary]]" << std::endl;
  fileOut << "summary = 'smoke'" << std::endl;
}

void
sansmic::Model::write_out_timestep_per_cell (int i, double p1, double p2)
{
  fileOut << "[[write_out_timestep_per_cell]]" << std::endl;
  fileOut << "i = " << i << std::endl;
  fileOut << "p1 = " << p1 << std::endl;
  fileOut << "p2 = " << p2 << std::endl;
}

void
sansmic::Model::write_out_timestep_totals (double qo)
{
  fileOut << "[[write_out_timestep_totals]]" << std::endl;
  fileOut << "qo = " << qo << std::endl;
}

void
sansmic::Model::write_out_timestep_insolubles (double p1, double p2)
{
  fileOut << "[[write_out_timestep_insolubles]]" << std::endl;
  fileOut << "p1 = " << p1 << std::endl;
  fileOut << "p2 = " << p2 << std::endl;
}

void
sansmic::Model::write_out_timestep_removed (double p1)
{
  fileOut << "[[write_out_timestep_removed]]" << std::endl;
  fileOut << "p1 = " << p1 << std::endl;
}

sansmic::JetModel::JetModel (void) { version = 1; }

sansmic::JetModel::JetModel (int ver) { version = ver; }

void
sansmic::JetModel::set_version (int ver)
{
  version = ver;
}

int
sansmic::JetModel::get_version (void)
{
  return version;
}

double
sansmic::JetModel::velocity (double Q, double r)
{
  return Q;
}

double
sansmic::JetModel::length (double u)
{
  return u;
}

sansmic::ODESolver::ODESolver (void) { neqn = 0; }

sansmic::ODESolver::ODESolver (int n, Derivable *f)
{
  neqn = n;
  dy_dx = f;
  yy = std::vector<double> (neqn + 1, 0.0);
  wt = std::vector<double> (neqn + 1, 0.0);
  p = std::vector<double> (neqn + 1, 0.0);
  yp = std::vector<double> (neqn + 1, 0.0);
  ypout = std::vector<double> (neqn + 1, 0.0);
  for (int i = 0; i < neqn + 1; i++)
    {
      std::vector<double> tmp = std::vector<double> (17, 0.0);
      phi.push_back (tmp);
    }
  alpha = std::vector<double> (13, 0.0);
  beta = std::vector<double> (13, 0.0);
  sig = std::vector<double> (14, 0.0);
  v = std::vector<double> (13, 0.0);
  w = std::vector<double> (13, 0.0);
  g = std::vector<double> (14, 0.0);
  psi = std::vector<double> (13, 0.0);
  iphase = 1;
  phase1 = true;
  x = 0;
  h = 0;
  hold = 0;
  istart = 1;
  start = true;
  told = 0;
  delsgn = 0;
  iwork1 = 0;
  iwork2 = 0;
  nornd = false;
  iwork3 = 0;
  iwork4 = 0;
  iwork5 = 0;
}

int
sansmic::ODESolver::num_eqn (void)
{
  return neqn;
}

void
sansmic::ODESolver::ode (std::vector<double> &y, double &t, double tout,
                         double &relerr, double &abserr, int &iflag)
{
  return;
}

void
sansmic::ODESolver::de (vector<double> &y, double &t, double tout,
                        double &relerr, double &abserr, int &iflag,
                        bool &phase1, bool &start, int &ns, bool &nornd,
                        int &k, int &kold, int &isnold)
{
  return;
}

void
sansmic::ODESolver::step1 (double &eps, bool &start, int &k, int &kold,
                           bool &crash, bool &phase1, int &ns, bool &nornd)
{
  return;
}

void
sansmic::ODESolver::intrp (double xout, vector<double> &yout, int kold)
{
  return;
}

sansmic::PlumeRise::PlumeRise (double delta_z, double alpha_entr,
                               std::vector<double> &conc)
{
  dz = delta_z;
  alpha = alpha_entr;
  co = conc;
  solver = ODESolver (neqn, (Derivable *)this);
}

void
sansmic::PlumeRise::func (double &x, std::vector<double> &y,
                          std::vector<double> &yp)
{
  yp[1] = 1.0;
  yp[2] = 2.0;
  yp[3] = 3.0;
}

sansmic::Salt::Salt (void)
{
  C_sat = 1.2019; // sg
  rho_s = 2.16;   // g/cm^3
  a_[0] = 1.0;
  a_[1] = 2.0;
  a_[2] = 3.0;
  a_[3] = 4.0;
  a_[4] = 5.0;
  a_[5] = 6.0;
  c_[0] = 1.0;
  c_[1] = 2.0;
  c_[2] = 3.0;
}

sansmic::Salt::Salt (double sg_max, double rho_solid)
{
  C_sat = sg_max;
  rho_s = rho_solid;
  a_[0] = 1.0;
  a_[1] = 2.0;
  a_[2] = 3.0;
  a_[3] = 4.0;
  a_[4] = 5.0;
  a_[5] = 6.0;
  c_[0] = 1.0;
  c_[1] = 2.0;
  c_[2] = 3.0;
}

void
sansmic::Salt::set_sg_max (double sg_max)
{
  C_sat = sg_max;
}

double
sansmic::Salt::get_sg_max (void)
{
  return C_sat;
}

void
sansmic::Salt::set_solid_density (double rho_solid)
{
  rho_s = rho_solid;
}

double
sansmic::Salt::get_solid_density (void)
{
  return rho_s;
}

void
sansmic::Salt::set_recession_rate_coeff (std::array<double, 6> new_coeff_a)
{
  for (int i = 0; i < 6; i++)
    {
      a_[i] = new_coeff_a[i];
    }
}

void
sansmic::Salt::set_density_conversion_coeff (std::array<double, 3> new_coeff_c)
{
  for (int i = 0; i < 3; i++)
    {
      c_[i] = new_coeff_c[i];
    }
}

std::array<double, 6>
sansmic::Salt::get_recession_rate_coeff ()
{
  return a_;
}

std::array<double, 3>
sansmic::Salt::get_sg_wt_pct_convert_coeff ()
{
  return c_;
}

double
sansmic::Salt::get_wt_pct (double sg, double temp)
{
  return sg;
}

double
sansmic::Salt::get_sg (double wt_pct, double temp)
{
  return wt_pct;
}

double
sansmic::Salt::get_recession_rate (double sg)
{
  return sg;
}

void
sansmic::trigad (int ns, int nf, std::vector<double> &y,
                 std::vector<double> &a, std::vector<double> &b,
                 std::vector<double> &c, std::vector<double> &d)
{
  return;
}

double
sansmic::sq (double val)
{
  return val * val;
}

double
sansmic::dabs (double d1)
{
  if (d1 >= 0.0)
    return d1;
  return -d1;
}

double
sansmic::ddim (double d1, double d2)
{
  if (d1 > d2)
    return d1 - d2;
  return 0.0;
}

double
sansmic::sign (double d1, double d2)
{
  if (d2 >= 0)
    return sansmic::dabs (d1);
  return -sansmic::dabs (d1);
}

int
sansmic::sign (int i1, int i2)
{
  if (i2 >= 0)
    return std::abs (i1);
  return -std::abs (i1);
}

double
sansmic::xnterp_dkr (double xf, std::vector<double> x, std::vector<double> y,
                     int nmax)
{
  double yf = -999.0;
  for (int i = 0; i < nmax; i++)
    {
      if (xf < x[i] && xf >= x[i + 1])
        {
          yf = y[i] + (xf - x[i]) / (x[i + 1] - x[i]) * (y[i + 1] - y[i]);
          break;
        }
    }
  return yf;
}

double
sansmic::bound (double value, double low, double high)
{
  return std::min (std::max (value, low), high);
}

sansmic::Model::Model (string out_prefix)
{
  salt = sansmic::Salt ();
  prefix = out_prefix;
  results = sansmic::Results ();
  jet_model = sansmic::JetModel ();
  stageNum = 0;
  n_cells = 0;
  n_nodes = 0;
  n_wells = 1;
  b_running = false;
  days = -1.0;
  timet = -1.0;
  V_tot = -1.0;
  f_dis0 = 1.0;
  stepNum = 0;
  C_cavAve = 0;
  Q_iTot = 0.0;
  Q_fTot = 0.0;
  b_injecting = false;
  injCell = 0;
  prodCell = 0;
  obiCell = 0;
  jetPlumeCell = 0;
  L_jet = 0;
  r_inj0 = 0;
  u_inj0 = 0;
  f_insol0 = 0.0;
  h_max = -1;
  floor_depth = 0;
  initialized = false;
  set_diffusion_beta_coeff (0.147);
  set_entrainment_coeff (0.09);
  set_molecular_diffusion_coeff (5.03e-5);
  set_eddy_coeff (1.142e5);
  set_plume_model_version (1);
  set_temperature_model_version (0);
  set_ODESolver_relative_tolerance (1.0e-4);
  set_ODESolver_absolute_tolerance (1.0e-2);
}

int
sansmic::Model::add_stage (Stage stage)
{
  stages.push_back (stage);
  return stages.size ();
}

int
sansmic::Model::get_num_cells (void)
{
  return n_cells;
}

void
sansmic::Model::set_num_cells (int n)
{
  if (initialized)
    throw std::exception ("Already initialized");
  n_cells = n;
  n_nodes = n + 1;
}

double
sansmic::Model::get_cavern_height (void)
{
  return h_max;
}

void
sansmic::Model::set_cavern_height (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  h_max = value;
}

double
sansmic::Model::get_initial_floor_depth (void)
{
  return floor_depth;
}

void
sansmic::Model::set_initial_floor_depth (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  floor_depth = value;
}

int
sansmic::Model::get_plume_model_version (void)
{
  return plumeModelVer;
}

void
sansmic::Model::set_plume_model_version (int ver)
{
  if (initialized)
    throw std::exception ("Already initialized");
  plumeModelVer = ver;
}

int
sansmic::Model::get_temperature_model_version (void)
{
  return tempModelVer;
}

void
sansmic::Model::set_temperature_model_version (int ver)
{
  if (initialized)
    throw std::exception ("Already initialized");
  tempModelVer = ver;
}

int
sansmic::Model::get_jet_model_version (void)
{
  return jet_model.get_version ();
}

void
sansmic::Model::set_jet_model_version (int ver)
{
  jet_model.set_version (ver);
  if (initialized)
    throw std::exception ("Already initialized");
}

double
sansmic::Model::get_entrainment_coeff (void)
{
  return alpha;
}

void
sansmic::Model::set_entrainment_coeff (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  alpha = value;
}

double
sansmic::Model::get_diffusion_beta_coeff (void)
{
  return beta;
}

void
sansmic::Model::set_diffusion_beta_coeff (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  beta = value;
}

double
sansmic::Model::get_molecular_diffusion_coeff (void)
{
  return D_mol;
}

void
sansmic::Model::set_molecular_diffusion_coeff (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  D_mol = value;
}

double
sansmic::Model::get_eddy_coeff (void)
{
  return D_0;
}

void
sansmic::Model::set_eddy_coeff (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  D_0 = value;
}

double
sansmic::Model::get_fraction_insolubles (void)
{
  return f_insol0;
}

void
sansmic::Model::set_fraction_insolubles (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  f_insol0 = value;
}

double
sansmic::Model::get_dissolution_factor (void)
{
  return f_dis0;
}

void
sansmic::Model::set_dissolution_factor (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  f_dis0 = value;
}

double
sansmic::Model::get_ODESolver_relative_tolerance (void)
{
  return relerr;
}

void
sansmic::Model::set_ODESolver_relative_tolerance (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  relerr = value;
}

double
sansmic::Model::get_ODESolver_absolute_tolerance (void)
{
  return abserr;
}

void
sansmic::Model::set_ODESolver_absolute_tolerance (double value)
{
  if (initialized)
    throw std::exception ("Already initialized");
  abserr = value;
}

void
sansmic::Model::init_model (void)
{
  sansmic::Model *self = this;
  x = 0;
  u = 0;
  initialized = true;

  temp = 75.0;

  C_hat = salt.get_sg_max ();
  C_wall = salt.get_solid_density ();
  w_hat = salt.get_wt_pct (C_hat, temp);

  b_running = false;
  b_injecting = true;
  b_times_up = false;
  b_obi_below_roof = true;
  b_no_rest_period = false;
  b_first_stage = true;

  idqi = 0;
  idqf = 1;
  nwrite1 = 0;
  nwrite2 = 0;
  dayOld = -1.0;
  stageNum = 0;

  V_insol = 0.0;
  volInsolVent = 0.0;
  h_insol = 0.0;
  err_cfac = 1.0;
  ppp = 0.777;
  Q_iOld = 0.0;
  Q_iTot = 0.0;
  Q_fOld = 0.0;
  Q_fTot = 0.0;
  Q_out = 0.0;
  p1 = 0;
  p2 = 0;
  h_insol = 0;
  izbs = 0;
  days = timet = 0.0;
  dz = 0.0;
  t_tot = 0.0;
  dt = 0.0;
  ntime = 0;
  diffCoeff = D_mol;
  Q_iTot = 0.0;
  Q_fTot = 0.0;
  timet = 0.0;
  t_last = 0;
  C_tmp = vector<double> (n_nodes + 1, 0.0);
  tunc = vector<double> (n_cells + 1, -1e10);
  aa = vector<double> (n_cells + 1, 0.0);
  voldkr = vector<double> (n_nodes + 1, 0.0);
  rcscr = vector<double> (n_nodes + 1, 0.0);
  r_cav0 = vector<double> (n_nodes + 1, 0.0);
  r_cav = vector<double> (n_nodes + 1, 0.0);
  h_cav = vector<double> (n_nodes + 1, 0.0);
  C_plume = vector<double> (n_nodes + 1, 0.0);
  u_plume = vector<double> (n_nodes + 1, 0.0);
  r_plume = vector<double> (n_nodes + 1, 0.0);
  z_cav = vector<double> (n_nodes + 1, 0.0);
  p2d = vector<double> (n_nodes + 1, 0.0);
  results.r_0 = vector<double> (n_nodes, 0.0);
  results.h_0 = vector<double> (n_nodes, 0.0);
  results.z_0 = vector<double> (n_nodes, 0.0);
  init_stage ();
  initialized = true;
  open_outfiles (false);
  write_tst_header ();
}

int
sansmic::Model::init_stage (void)
{
  return stageNum;
}

int
sansmic::Model::end_stage (void)
{
  write_log_end_stage ();
  return 1;
}

void
sansmic::Model::run_sim (void)
{
  init_model ();
  close_outfiles ();
}

int
sansmic::Model::run_stage ()
{
  return 1;
}

int
sansmic::Model::run_step ()
{
  return 0;
}

int
sansmic::Model::leach ()
{
  return stepNum;
}

sansmic::Results
sansmic::Model::get_current_state (void)
{
  sansmic::Results new_results = sansmic::Results ();
  new_results.r_0 = vector<double> (n_nodes, 0.0);
  new_results.h_0 = vector<double> (n_nodes, 0.0);
  new_results.z_0 = vector<double> (n_nodes, 0.0);
  for (int j = 0; j < n_nodes; j++)
    {
      new_results.r_0[j] = results.r_0[j];
      new_results.h_0[j] = results.h_0[j];
      new_results.z_0[j] = results.z_0[j];
    }

  save_results (new_results, false);
  return new_results;
}

void
sansmic::Model::save_results (sansmic::Results &new_results, bool to_file)
{
  new_results.t.push_back (t_tot);
  new_results.dt.push_back (dt);
  new_results.step.push_back (stepNum);
  new_results.V_cavTot.push_back (V_tot);
  new_results.err.push_back (err_cfac);
  new_results.sg_cavAve.push_back (C_cavAve);
  new_results.V_injTot.push_back (Q_iTot);
  new_results.V_fillTot.push_back (Q_fTot);
  new_results.stage.push_back (stageNum);
  new_results.phase.push_back ((int)b_injecting);
  new_results.injCell.push_back (injCell);
  new_results.prodCell.push_back (prodCell);
  new_results.obiCell.push_back (obiCell);
  new_results.plmCell.push_back (jetPlumeCell);
  new_results.z_plm.push_back (z_cav[1] - h_cav[jetPlumeCell]);
  new_results.z_inj.push_back (z_cav[1] - h_cav[injCell]);
  new_results.z_prod.push_back (z_cav[1] - h_cav[prodCell]);
  new_results.l_jet.push_back (L_jet);
  new_results.r_jet.push_back (r_inj0);
  new_results.u_jet.push_back (u_inj0);

  // DO OUTPUTS
  if (to_file)
    {
      write_out_timestep_summary ();
    }

  for (int i = 1, j = 1; j <= n_nodes; j++)
    {
      p1 = i + 1;
      p2 = j + 2;
      if (to_file)
        {
          write_out_timestep_per_cell (i, p1, p2);
        }
    }
  vector<double> _r_cav, _dr_cav, _sg, _theta, _Q_inj, _V, _f_dis, _f_flag,
      _xincl, _amd, _D_coeff, _dC_dz, _C_old, _C_new, _dC_dt, _dr_dt, _C_plm,
      _u_plm, _r_plm;
  for (int i = 1; i <= n_nodes; i++)
    {
      p1 = _pi_ * dz * sq (r_cav[i]) * cf_to_bbl;
      p2 = V_injSigned[i] * cf_to_bbl * per_h_to_per_d;
      _r_cav.push_back (r_cav[i]);
      _dr_cav.push_back (r_cav[i] - r_cav0[i]);
      _sg.push_back (C_cav[i]);
      _theta.push_back (phi[i]);
      _Q_inj.push_back (p2);
      _V.push_back (p1);
      _f_dis.push_back (f_dis_prt[i]);
      _f_flag.push_back (f_disType[i]);
      _xincl.push_back (x_incl[i]);
      _amd.push_back (amd_prt[i] * cf_to_bbl * per_h_to_per_d);
      _D_coeff.push_back (akd_prt[i]);
      _dC_dz.push_back (ca_prt[i]);
      _C_old.push_back (C_tmp[i]);
      _C_new.push_back (C_cav[i]);
      _dC_dt.push_back (dC[i]);
      _dr_dt.push_back (dr_prt[i]);
      _C_plm.push_back (C_plume[i]);
      _u_plm.push_back (u_plume[i]);
      _r_plm.push_back (r_plume[i]);
    }
  new_results.r_cav.push_back (_r_cav);
  new_results.dr_cav.push_back (_dr_cav);
  new_results.sg.push_back (_sg);
  new_results.theta.push_back (_theta);
  new_results.Q_inj.push_back (_Q_inj);
  new_results.V.push_back (_V);
  new_results.f_dis.push_back (_f_dis);
  new_results.f_flag.push_back (_f_flag);
  new_results.xincl.push_back (_xincl);
  new_results.amd.push_back (_amd);
  new_results.D_coeff.push_back (_D_coeff);
  new_results.dC_dz.push_back (_dC_dz);
  new_results.C_old.push_back (_C_old);
  new_results.C_new.push_back (_C_new);
  new_results.dC_dt.push_back (_dC_dt);
  new_results.dr_dt.push_back (_dr_dt);
  new_results.C_plm.push_back (_C_plm);
  new_results.u_plm.push_back (_u_plm);
  new_results.r_plm.push_back (_r_plm);

  if (to_file)
    {
      write_out_timestep_totals (Q_outBPD);
    }
  new_results.Q_out.push_back (Q_outBPD);
  new_results.sg_out.push_back (C_cav[prodCell]);

  if (to_file)
    {
      write_out_timestep_insolubles (p1, p2);
    }
  new_results.V_insolTot.push_back (p1);
  new_results.V_insolVent.push_back (p2);
  new_results.h_insol.push_back (h_insol);
  new_results.z_insol.push_back (z_cav[1] - h_insol);
  new_results.z_obi.push_back (z_cav[1] - h_bkt);

  if (to_file)
    {
      write_out_timestep_removed (p1);
    }
}

double
sansmic::Model::brine_mass (double dz)
{
  return dz;
}

void
sansmic::Model::remove_salt (int i, double recession_rate, double dt,
                             double dz)
{
  return;
}

void
sansmic::Model::plume (double ci, double zi, double zb, double &x, double &u,
                       double dz, double alpha, int nmax, double &r,
                       double &cpl)
{
  return;
}

void
sansmic::Model::slope (int i)
{
  return;
}
