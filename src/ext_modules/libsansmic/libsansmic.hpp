// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file libsansmic.hpp
 * @brief libsansmic API definitions
 */

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#ifndef SANSMIC_H
#define SANSMIC_H

using namespace std;

namespace sansmic {

// physical constants and conversion factors
const double _pi_ = 3.141592653589793;        //!< pi without having to import
const double std_gravity = 9.80665;           //!< standard gravity (m/s²)
const double inch = 0.0254;                   //!< 1 in. := *y* m
const double foot = 0.3048;                   //!< 1 ft := *y* m
const double in_to_ft = 1 / 12.0;             //!< ft/in
const double ft_to_in = 12.0;                 //!< in/ft
const double square_inch = 0.00064516;        //!< 1 in² := *y* m²
const double square_foot = 0.09290304;        //!< 1 ft² := *y* m²
const double cubic_inch = 1.6387064e-05;      //!< 1 in³ := *y* m³
const double liquid_gallon = 0.003785411784;  //!< 1 gal := *y* m³
const double cubic_foot = 0.028316846592;     //!< 1 ft³ := *y* m³
const double oil_barrel = 0.158987294928;     //!< 1 bbl := *y* m³
const double bbl_to_cf = 539.0 / 96.0;        //!< 1 bbl := *y* ft³
const double cf_to_bbl = 96.0 / 539.0;        //!< 1 ft³ := *y* bbl
const double pound_mass = 0.45359237;         //!< 1 lb := *y* kg
const double pound_force = 4.4482216152605;   //!< 1 lbf := *y* N
const double degree_Fahrenheit = 5.0 / 9.0;   //!< 1 F° := *y* K
const double zero_Fahrenheit = 255.15 + 2.0 / 9.0;  //!< 0 °F := *y* K
const double minute = 60.0;                         //!< s / min
const double hour = 3600.0;                         //!< s / h
const double day = 86400.0;                         //!< s / d
const double julian_year = 86400.0 * 365.25;        //!< s / aⱼ
const double per_d_to_per_h = hour / day;           //!< h / d
const double per_h_to_per_d = day / hour;           //!< d / h

// error codes
const int MISSING_DAT_FILE = 55;  //!< unable to open the "dat" input file
const int INVALID_LOG_FILE = 54;  //!< unable to open the "log" log file
const int INVALID_OUT_FILE = 56;  //!< unable to open the "out" output file
const int INVALID_TST_FILE = 59;  //!< unable to open the "tst" file
const int NDIV_CHANGES_THROUGHOUT_FILE = 101;  //!< inconsistent n_div value
const int OBI_AT_CAVERN_FLOOR = 201;  //!< the OBI has reached the cavern floor
const int ODE_IFLAG_SIX = 306;        //!< ODE solver failed with value 6
const int ODE_IFLAG_SEVEN = 307;      //!< ODE solver failed with value 7
const int ODE_BAD_ABSERR = 300;  //!< ODE solver failed due to bad error value
const int UNIMPLEMENTED_FLOW_TABLES =
    900;  //!< flow tables not currently implemented
const int UNIMPLEMENTED_GEOMETRY_IDATA =
    901;  //!< the geometry type is not implemented yet

// ENUM definitions

/// @brief The leaching operations mode for sansmic.
enum LeachMode : int {
  Ordinary = 0,    //!< fill with raw water, produce brine
  Withdrawal = 1,  //!< fill with raw water, produce oil
  LeachFill = 2,   //!< fill with oil and raw water, produce brine
  OilFill = -1     //!< fill with oil, produce brine
};

/// @brief The format of the data; will go away once python code fully handles
/// these.
enum GeomFormat : int {
  RadiusList = 0,    //!< equally spaced radius measurements, from bottom to top
  VolumeList = 1,    //!< depth data followed by corresponding partial volumes
  VolumeTable = -1,  //!< list of depth-volume pairs
  RadiusTable = 2    //!< list of depth-radius pairs
};

// utility mathematics functions in "util.cpp"
double sq(double val);
double dabs(double d1);
double ddim(double d1, double d2);
double sign(double d1, double d2);
int sign(int i1, int i2);
double xnterp_dkr(double xf, vector<double> x, vector<double> y, int nmax);
double bound(double value, double low, double high);
void trigad(int ns, int nf, vector<double> &y, vector<double> &a,
            vector<double> &b, vector<double> &c, vector<double> &d);

// Stage struct definition in "stage.cpp"
/// @brief Instructions and input values for a single stage of SANSMIC.
struct Stage {
  Stage();
  void debug_log(ofstream &fout, int stageNum = 0);

  string title;  //!< stage title

  LeachMode leach_mode;  //!< the run mode (leaching mode) to use

  int print_freq;     //!< the number of timesteps (*not hours*) between output
                      //!< steps
  int is_subsequent;  //!< false if this is the fisrt stage, true if a
                      //!< subsequent stage

  double t_rest;  //!< the number of hours to wait at the end of the stage
  double stop_cond_val;  //!< volume (>0, in bbl) or OBI (<0, in ft) to use as
                         //!< end-of-stage, 0 for time based
  double h_inj;    //!< height of the injection point above the floor (in feet)
  double h_prod;   //!< height of the production point above the floor (in feet)
  double h_obi;    //!< height of the OBI above the floor (in feet)
  double Q_raw;    //!< injcetion rate (in bbl/hr)
  double r_tbgID;  //!< radius (in inches) of the inside of the inner string
  double r_tbgOD;  //!< radius (in inches) of the outside of the inner string
  double r_csgID;  //!< radius (in inches) of the inside of the outer string
  double r_csgOD;  //!< radius (in inches) of the outside of the outer string
  double sg_raw;   //!< specific gravity of injected water
  double sg_init;  //!< specific gravity of the cavern at the start of the stage
  double dt;       //!< timestep (in hours)
  double t_stage;  //!< injection duration (in hours)
  double Q_oil;    //!< fill rate for LeachFill mode (in bbl/hr)
};

// Scenario struct definitions in "scenario.cpp"
/// @brief Scenario options that are constant through all stages.
struct Scenario {
  Scenario();

  int add_stage(Stage stage);

  void debug_log(ofstream &fout);

  string title;     //!< a title for the full scenario
  string comments;  //!< longer comments about scenario

  GeomFormat geometry_format;  //!< geometry input format

  int num_cells;  //!< the number of model cells

  double cavern_height;    //!< the height of the cavern (dz = height/num_cells)
  double floor_depth;      //!< the depth (MD) of the initial cavern TD
  double ullage_standoff;  //!< the height above TD used to calculate ullage
  double fraction_insolubles;  //!< the fraction of insolubles

  vector<double> geom_radii;    //!< list of radii (in feet)
  vector<double> geom_depths;   //!< list of depths (in feet; if applicable)
  vector<double> geom_volumes;  //!< list of volumes (in bbl; if applicable)

  vector<Stage> stages;  //!< simulation stages

  // advanced settings;
  int jet_model_version;          //!< the ver for the jet model (0 = off)
  int plume_model_version;        //!< the ver for the plume model (always 1)
  int temperature_model_version;  //!< the ver for the temp model (always 0)

  double max_brine_sg;         //!< maximum brine density
  double solid_density;        //!< solid density of salt rock
  double entrainment_coeff;    //!< the alpha entrainment coefficient
  double diffusion_beta;       //!< the beta diffusion coefficient
  double molecular_diffusion;  //!< the D_mol diffusion coefficient
  double eddy_coefficient;     //!< the diffusion eddy coefficient D_0
  double relative_error;       //!< the relative tolerance for the ODE solver
  double absolute_error;       //!< the absolute tolerance for the ODE solver

  // deprecated options
  int coallescing_wells;      //!< the number of identical wells coallescing
  double well_separation;     //!< the separation between coallescing wells
  double dissolution_factor;  //!< the dissolution factor (this should be 1.0!)
};

// Salt definition in "salt.cpp"
/// @brief Class defining the properties for a certain type of rock salt.
class Salt {
 public:
  Salt(void);
  Salt(double sg_max, double rho_solid);

  double get_sg_max(void);
  double get_solid_density(void);
  array<double, 6> get_recession_rate_coeff();
  array<double, 3> get_sg_wt_pct_convert_coeff();

  void set_sg_max(double sg_max);
  void set_solid_density(double rho_solid);
  void set_recession_rate_coeff(array<double, 6> new_coeff_a);
  void set_density_conversion_coeff(array<double, 3> new_coeff_c);

  double wt_pct(double sg, double temp = 75.0);
  double sg(double wt_pct, double temp = 75.0);
  double recession_rate(double sg);

 private:
  double C_sat;  //!< specific gravity of saturated brine
  double rho_s;  //!< density of the surrounding rock in g/cm^3

  array<double, 6> a_;  //!< the coeff for wall recession rate formula
  array<double, 3> c_;  //!< the coeff for converting b/w sg & wt pct
};

// JetModel definition in "jetmodel.cpp"
/// @brief Model of the water jet flow out the end of tubing.
class JetModel {
 public:
  JetModel(void);
  JetModel(int ver);

  void set_version(int ver);
  int get_version(void);

  double velocity(double Q, double r);
  double length(double u);

 private:
  int version;  //!< the jet model version - 0 is off
};

// Interface definition local
/// @brief Interface for a class that provides n-derivatives for an ODE func.
class Derivable {
 public:
  /**
   * @brief A function that provides the necessary derivatives for the ODE
   * solver.
   * @param x the value to calculate derivatives at
   * @param y the f(y) values
   * @param yp the f'(y) values
   */
  virtual void func(double &x, vector<double> &y, vector<double> &yp) = 0;
  int neqn;  //!< The number of equations being calculated
};

// ODESolver definition in "odesolver.cpp"
/**
 * @brief The ODE solver class
 * @details The function of this code is completely explained and documented in
 * the text, *Computer solution of ordinary differential equations, the initial
 * value problem* by L. F. Shampine and M. K. Gordon. Further details on use of
 * this code are available in *Solving ordinary differential equations with
 * ODE, STEP, and INTRP*, by L. F. Shampine and M. K. Gordon, SLA-73-1060.
 */
class ODESolver {
 public:
  ODESolver(void);
  ODESolver(int n, Derivable *f);
  int num_eqn(void);
  void ode(vector<double> &y, double &t, double tout, double &relerr,
           double &abserr, int &iflag);

 protected:
  void de(vector<double> &y, double &t, double tout, double &relerr,
          double &abserr, int &iflag, bool &phase1, bool &start, int &ns,
          bool &nornd, int &k, int &kold, int &isnold);
  void intrp(double xout, vector<double> &yout, int kold);
  void step1(double &eps, bool &start, int &k, int &kold, bool &crash,
             bool &phase1, int &ns, bool &nornd);

 protected:
  int neqn;  //!< the number of equations

  Derivable *dy_dx;  //!< the derivative function object

 private:
  bool phase1;  //!< internal work variable
  bool start;   //!< internal work variable
  bool nornd;   //!< internal work variable

  int iphase;  //!< internal work variable
  int istart;  //!< internal work variable
  int iwork1;  //!< internal work variable
  int iwork2;  //!< internal work variable
  int iwork3;  //!< internal work variable
  int iwork4;  //!< internal work variable
  int iwork5;  //!< internal work variable

  double x;       //!< internal work variable
  double h;       //!< internal work variable
  double hold;    //!< internal work variable
  double told;    //!< internal work variable
  double delsgn;  //!< internal work variable

  vector<double> yy;     //!< internal work variable
  vector<double> wt;     //!< internal work variable
  vector<double> p;      //!< internal work variable
  vector<double> yp;     //!< internal work variable
  vector<double> ypout;  //!< internal work variable
  vector<double> alpha;  //!< internal work variable
  vector<double> beta;   //!< internal work variable
  vector<double> sig;    //!< internal work variable
  vector<double> v;      //!< internal work variable
  vector<double> w;      //!< internal work variable
  vector<double> g;      //!< internal work variable
  vector<double> psi;    //!< internal work variable

  vector<vector<double> > phi;  //!< internal work variable
};

// PlumeRise definition in "plumerise.cpp"
/// @brief Class defining the function for the plume rise integration in ODE.
class PlumeRise : public Derivable {
 public:
  PlumeRise(double delta_z, double alpha_entr, vector<double> &conc);

  void func(double &x, vector<double> &y, vector<double> &yp);

  ODESolver solver;    //!< solver for the ode
  const int neqn = 3;  //!< number of equations used

 private:
  double dz;     //!< cell size in ft
  double alpha;  //!< entrainment coefficient

  vector<double> co;  //!< concentration vector for the cavern
};

// Results struct definition is local
/// @brief Structure containing the results of the simulation
struct Results {
  /* vectors by height */
  vector<double> r_0;  //!< initial radius by node
  vector<double> z_0;  //!< node depth
  vector<double> h_0;  //!< node height

  /* vectors by time */
  vector<int> step;      //!< the step number
  vector<int> stage;     //!< stage number
  vector<int> phase;     //!< phase 1 (inject), phase 0 (static)
  vector<int> injCell;   //!< cell containing EOT
  vector<int> prodCell;  //!< cell containing production
  vector<int> obiCell;   //!< cell containing OBI
  vector<int> plmCell;   //!< cell containing plume stagnation

  vector<double> t;            //!< time of output
  vector<double> err;          //!< mass ballance ratio
  vector<double> z_obi;        //!< obi _depth_
  vector<double> z_inj;        //!< injection point _depth_
  vector<double> z_prod;       //!< production point _depth_
  vector<double> z_plm;        //!< plume stagnation _depth_
  vector<double> z_insol;      //!< insoluble top _depth_
  vector<double> h_insol;      //!< insoluble _height_
  vector<double> l_jet;        //!< injection jet length
  vector<double> r_jet;        //!< injection radius
  vector<double> u_jet;        //!< injection point velocity
  vector<double> V_injTot;     //!< total injected volume
  vector<double> V_fillTot;    //!< total fill/withdraw (oil) volume
  vector<double> V_cavTot;     //!< total cavern volume
  vector<double> V_insolTot;   //!< total insolubles created
  vector<double> V_insolVent;  //!< volume insolubles vented
  vector<double> Q_out;        //!< current production rate
  vector<double> sg_out;       //!< outlet specific gravity
  vector<double> sg_cavAve;    //!< cavern average specific gravity
  vector<double> dt;           //!< timestep (in case stage changes)

  /* data by time and cell */
  vector<vector<double> > r_cav;    //!< cavern radius
  vector<vector<double> > dr_cav;   //!< change in radius
  vector<vector<double> > sg;       //!< specific gravity
  vector<vector<double> > theta;    //!< wall angle
  vector<vector<double> > Q_inj;    //!< injection volume
  vector<vector<double> > V;        //!< cell volume
  vector<vector<double> > f_dis;    //!< dissolution factor
  vector<vector<double> > f_flag;   //!< factor type flag
  vector<vector<double> > xincl;    //!< wall angle correction
  vector<vector<double> > amd;      //!< debug
  vector<vector<double> > D_coeff;  //!< diffusion coefficient
  vector<vector<double> > dC_dz;    //!< change in sg vertically
  vector<vector<double> > C_old;    //!< previous SG
  vector<vector<double> > C_new;    //!< current SG
  vector<vector<double> > dC_dt;    //!< C_old - C_new
  vector<vector<double> > dr_dt;    //!< recession rate
  vector<vector<double> > C_plm;    //!< plume sg
  vector<vector<double> > u_plm;    //!< plume veloctiy
  vector<vector<double> > r_plm;    //!< plume radius
};

// The pieces of the model that are accessed by multiple functions
/// @brief A base model class
class BaseModel {
 public:
  bool get_is_running(void);
  int get_num_stages(void);
  int get_current_stage(void);
  double get_current_time(void);
  double get_current_volume(void);
  vector<Stage> get_stages(void);
  void set_verbosity_level(int verb);
  void set_use_tstfile(bool use_file);
  void set_use_outfile(bool use_file);
  Results get_results(void);
  Results get_current_state(void);

 protected:
  string prefix;         //!< output file prefix
  Results results;       //!< results object
  vector<Stage> stages;  //!< simulation stages

  Salt *salt = NULL;             //!< salt properties object pointer
  JetModel *jet_model = NULL;    //!< the jet model object pointer
  PlumeRise *plume_rise = NULL;  //!< plume rise model object pointer

  bool b_is_injecting = true;   //!< is the stage injecting or static?
  bool b_is_running = false;    //!< continue running the timestep?
  bool b_use_outfile = false;   //!< create and use an .OUT file?
  bool b_use_tstfile = false;   //!< create and use a .TST file?
  double abserr = 1.0e-2;       //!< ODE solver: absolute tolerance
  double C_cavAve = 0.0;        //!< average cavern brine sg
  double C_inj;                 //!< specific gravity of injected water
  double days = 0.0;            //!< time: current loop time, in days
  double dt = 0.1;              //!< timestep size (h)
  double dz = 0.0;              //!< height of each cell
  double err_cfac = 1.0;        //!< convergence level
  double h_inj = 0.0;           //!< height: injection string EOT
  double h_insol = 0.0;         //!< height: of top of insolubles
  double h_obi = 0.0;           //!< height: initial OBI/blanket level
  double h_prd = 0.0;           //!< height: production string EOT
  double h_max = 0.0;           //!< height: total cavern height
  double L_jet = 0.0;           //!< jet model: jet length
  double Q_fOld = 0.0;          //!< oil fill: last timestep
  double Q_fTot = 0.0;          //!< oil fill: total volume
  double Q_iOld = 0.0;          //!< injection: last timestep
  double Q_iTot = 0.0;          //!< injection: total volume
  double Q_out = 0.0;           //!< out flow: rate
  double Q_outBPD = 0.0;        //!< out flow: rate in BPD
  double r_inj0 = 0.0;          //!< jet model: injection point radius
  double relerr = 1.0e-4;       //!< ODE solver: relative tolerance
  double slctim = 0.0;          //!<
  double t_tot = 0.0;           //!< time: total time, in hours
  double timet = 0.0;           //!< time: total time, in days
  double u_inj0 = 0.0;          //!< jet model: injection point velocity
  double V_insol = 0.0;         //!< volume: insolubles in cavern
  double V_insolVent = 0.0;     //!< volume: insolubles vented out
  double V_tot = 0.0;           //!< volume: total volume
  double V_ullage = 0.0;        //!< volume: ullage available
  double V_usable = 0.0;        //!< volume: total usable
  double volRemoved = 0.0;      //!< volume: volume salt removed
  int injCell = 0;              //!< cell containing the injection EOT
  int izbs = 0;                 //!< current OBI cell index
  int jetPlumeCell = 0;         //!< cell containing top of plume
  int n_nodes = 0;              //!< number of nodes
  int obiCell = 0;              //!< cell containing the interface
  int prodCell = 0;             //!< production string EOT cell
  int stageNum = 0;             //!< the current stage number
  int stepNum = 0;              //!< current step number
  int verbosity = 0;            //!< verbosity setting for output
  vector<double> akd_prt;       //!< AKD for printing
  vector<double> amd_prt;       //!< AMD for printing
  vector<double> C_cav;         //!< cavern brine sg vector
  vector<double> C_plume;       //!< plume concentration vector
  vector<double> C_tmp;         //!< concentration
  vector<double> ca_prt;        //!< CA for printing
  vector<double> dC;            //!< change in concentration between steps
  vector<double> dr_prt;        //!< delta radius (for output)
  vector<double> f_dis_prt;     //!< dissolution factor (for output)
  vector<double> h_cav;         //!< cell elevations from floor
  vector<double> phi;           //!< wall angle of the cell
  vector<double> r_cav;         //!< current cavern radius
  vector<double> r_cav0;        //!< initial cavern radius
  vector<double> r_plume;       //!< plume radius vector
  vector<double> rcscr;         //!<
  vector<double> tanTheta;      //!< the tangent of the wall angle
  vector<double> u_plume;       //!< plume velocity vector
  vector<double> V_injSigned;   //!<
  vector<double> V_saltRemove;  //!< volume salt removed
  vector<double> x_incl;        //!< wall angle factor
  vector<double> z_cav;         //!< cell measured depths
  vector<int> f_disType;        //!< dissolution regime indicator

  void save_results(Results &new_results);

  void write_version_info(ofstream &sout);
  void write_daily_column_header(ofstream &sout);
  void write_daily_summary(ofstream &sout, int stage, bool inject);
  void write_daily_end_of_phase(ofstream &sout, const char *text);
  void write_simulation_init(ofstream &sout);
  void write_stage_init(ofstream &sout);
  void write_stage_summary(ofstream &sout);
  void write_simulation_summary(ofstream &sout);
  void write_detailed_results(ofstream &sout);
};

// The sansmic numerical model is located in "model.cpp"
/// @brief The main SANSMIC model class.
class Model : public BaseModel {
 public:
  Model();
  Model(string out_prefix);
  void configure(Scenario scen);

  void run_sim(void);
  int run_stage();
  int run_step();
  int init_stage(void);
  int end_stage();
  int open_outfiles(bool append = false);
  void close_outfiles(void);

 private:
  void init_vars(void);

  string Q_fill_table_name;  //!< name of fill table file
  string Q_inj_table_name;   //!< name of injection table file

  ofstream fileOut;  //!< output file
  ofstream fileLog;  //!< log file
  ofstream fileTst;  //!< debug ".tst" file

  GeomFormat dataFmt;  //!< the format of the geometry data

  double alpha;    // entrainment coefficient
  double cpl;      //
  double r;        // radius
  double u;        // velocity
  double x;        //
  double m_brine;  //

  int nEqn;  //

  array<int, 2> idxq;   // USED FOR QI/QF TABLES
  array<int, 2> iqnum;  // USED FOR QI/QF TABLES

  bool b_first_stage;     // is this the first stage or a subsequent stage
  bool b_initialized;     // has the model been initialized
  bool b_just_saved;      // did we just save? don't save twice!
  bool b_obi_below_roof;  // did the obi exit the cavern
  bool b_times_up;        // has this step reached the end of the stage
  bool b_use_fill_table;  // use a flow table for product fill
  bool b_use_inj_table;   // use a flow table for brine injection

  double beta;          // diffusion beta coefficient
  double C_bar;         // average specific gravity in plume
  double C_cav0;        // initial specific gravity in cavern
  double C_hat;         // specific density saturated brine
  double C_plm;         // plume concentraion
  double C_wall;        // solid density
  double c1;            //
  double c2;            //
  double c3;            //
  double ca;            // temp variable A
  double cb;            // temp variable B
  double cc;            // temp variable C
  double cnext;         // next concentration
  double cold;          // old concentration
  double cpln;          //
  double D_0;           // eddy diffusion coefficient
  double D_mol;         // molecular diffusion coefficient
  double dayOld;        // used for determining daily output
  double dC_dt;         // concentration derivative w.r.t. time
  double diffCoeff;     //
  double dr_dt;         // recession rate
  double dt_dz2;        // dt / dz^2
  double dt_min;        // minimum recommended dt based on injection rates
  double dum;           // temp variable
  double duml;          // temp variable
  double dx_sep;        // separation between coallescing wells
  double dz_inc;        //
  double f_dis0;        // base adjustment of dissolution rate from input
  double f_insol0;      // fraction insolubles
  double fallf;         // fraction of falling insolubles
  double h_obiSave;     // OBI saved between stages
  double h_uso;         //
  double m_brineNew;    //
  double m_brineOld;    //
  double m_plume;       //
  double m_saltRemove;  //
  double p1, p2, p3, p4, p5, p7, p8, p9, p10;
  double Q_fill;         //
  double Q_fillBPD;      //
  double Q_fillSav;      //
  double Q_in;           //
  double Q_inBPD;        //
  double Q_iSav;         //
  double r_cavMax;       // maximum cavern diameter
  double r_cavMin;       // minimum cavern diameter
  double r_csgOD;        // outer radius of outer casing
  double r_innerPipeOD;  // outer radius of inner pipe
  double r_outerPipeID;  // inner radius of outer casing
  double r_sqd;          // um, volume not radius squared. Actually = r^2 dz
  double r_stag;         // radius of plume at plume stagnation level
  double r_tbgID;        // inner radius of inner pipe
  double rho_sat;        // saturated density
  double runMode;        // run mode
  double S_d;            //
  double sig;            //
  double stopCriteria;   //
  double t;              // current time
  double t_end;          // time to end injection
  double t_last;         //
  double t_wait;         // workover period duration
  double temp;           // temperature
  double theta;          // wall angle
  double V_cell;         //
  double V_insolRemain;  // insolubles within an cell
  double V_plume;        //
  double V_stop;         // stop when oil remaining less than value
  double V_used;         //
  double vdif;           //
  double vel;            //
  double velfac;         //
  double vflo;           // flow into plume
  double vi;             //
  double vip;            //
  double vppl;           // plume volume
  double vrn;            //
  double vtpl;           // total plume volume height
  double w;              //
  double w_d;            // downwind differencing variable
  double w_hat;          // wt pct. of saturated brine
  double w_u;            // upwind differencing variable
  double watsal;         //
  double watsol;         //
  double z_obi;          // obi depth
  double z_obi_stop;     // OBI-final (for stage termination criteria)
  double z_TD0;          //!< bottom of the cavern
  double z_ullage;       //

  int i_obi;              // temp blanket cell index
  int i_obiOld;           // temp blanket cell index
  int idqf;               //
  int idqi;               //
  int im;                 // i minus 1
  int injCellBelow;       //
  int ip;                 // i plus 1
  int jetPlumeCellBelow;  // index cell below jet
  int jplp;               // bounded jetPlumeCell + 1
  int jpls;               // jet plume cell
  int L;                  // TODO: ew ew ew change this var name
  int m;                  //
  int maxProdOrJet;       // highest cell that is either production or jet
  int minProdOrJet;       // lowest cell that is either production or jet
  int n_cells;            // number of cells
  int n_wells;            // number of coallescing caverns
  int noFlowCell;         //
  int obiCellBelow;       //
  int plumeModelVer;      // the version of the plume model
  int print_freq;         //
  int startCell;          //
  int tempModelVer;       // the version of the temperature model

  vector<double> aa;        //
  vector<double> cb_prt;    // used for printing CB
  vector<double> cosTheta;  //
  vector<double> depdkr;    // temp depth vector
  vector<double> f_dis;     // adjusted ratio of dissolution rate
  vector<double> f_disSav;  // base ratio of actual-to-assumed dissolution rate
  vector<double> f_insol;   // volume ratio of insolubles to salt
  vector<double> p2d;
  vector<double> qtim;
  vector<double> qval;
  vector<double> ss;
  vector<double> tunc;
  vector<double> vec_A;   // A-overline in (2.8)
  vector<double> vec_B;   // B-overline in (2.8)
  vector<double> vec_C;   // C-overline in (2.8)
  vector<double> vec_D;   // D-overline in (2.8)
  vector<double> voldkr;  // temp variable for ullage

  int leach();
  double brine_mass(double dz);
  void remove_salt(int i, double recession_rate, double dt, double dz);
  void plume(double ci, double zi, double zb, double &x, double &u, double dz,
             double alpha, int nmax, double &r, double &cpl);
  void slope(int i);
};

}  // namespace sansmic
#endif
