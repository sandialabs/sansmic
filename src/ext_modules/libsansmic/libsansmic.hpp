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

#ifndef VERSION_INFO
#define VERSION_INFO "local"
#endif

#ifndef SANSMIC_H
#define SANSMIC_H

using namespace std;

namespace sansmic {

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

/**
 * @brief Shortcut for square
 * @details Because FORTRAN and Python have an operator for power, a shortcut
 * function was written for brevity and to ensure there are no mistakes due to
 * typos between two complex groupings of products/sums that are then squared.
 * @param val
 * @return double
 */
double sq(double val);

/**
 * @brief absolute value for doubles
 * @param d1 value
 * @return double
 */
double dabs(double d1);

/**
 * @brief implementation of FORTRAN DDIM1 function in C++
 * @param d1 value 1
 * @param d2 value 2
 * @return double
 */
double ddim(double d1, double d2);

/**
 * @brief implementation of FORTRAN DSIGN function in C++
 * @param d1
 * @param d2
 * @return double
 */
double sign(double d1, double d2);

/**
 * @brief implementation of FORTRAN ISIGN function in C++
 * @param i1
 * @param i2
 * @return int
 */
int sign(int i1, int i2);

/**
 * @brief interpolation function
 * @param xf x value at which to interpolate
 * @param x known x values
 * @param y known y values
 * @param nmax length of arrays
 * @return double, interpolated y value at xf
 */
double xnterp_dkr(double xf, vector<double> x, vector<double> y, int nmax);

/**
 * @brief bound a value within a certain range
 * @param value original value
 * @param low bound on low end
 * @param high bound on high end
 * @return double bounded value
 */
double bound(double value, double low, double high);

/**
 * @brief Tridiagonal solver
 * @param ns lower boundary
 * @param nf last unknown point
 * @param y concentration (C_new in sansmic::Model)
 * @param a main-1 diagonal
 * @param b main diagonal
 * @param c main+1 diagonal
 * @param d right hand side
 */
void trigad(int ns, int nf, vector<double> &y, vector<double> &a,
            vector<double> &b, vector<double> &c, vector<double> &d);

/**
 * @brief The leaching/operations mode for SANSMIC.
 */
enum LeachMode : int {
  Ordinary = 0,    //!< fill with raw water, produce brine
  Withdrawal = 1,  //!< fill with raw water, produce oil
  LeachFill = 2,   //!< fill with oil and raw water, produce brine
  OilFill = -1     //!< fill with oil, produce brine
};

/**
 * @brief The format of the input data.
 */
enum GeomFormat : int {
  RadiusList = 0,    //!< equally spaced radius measurements, from bottom to top
  VolumeList = 1,    //!< depth data followed by corresponding partial volumes
  VolumeTable = -1,  //!< list of depth-volume pairs
  RadiusTable = 2    //!< list of depth-radius pairs
};

/**
 * @brief Instructions and input values for a single stage of SANSMIC.
 */
struct Stage {
  /**
   * @brief Stage constructor with default values.
   */
  Stage();

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

  /**
   * @brief Output debug data to a stream.
   * @param fout Open stream to output to
   */
  void debug_log(ofstream &fout);
};

/**
 * @brief Scenario options that are constant through all stages.
 */
struct Scenario {
  /**
   * Blank options with default version and coefficient values.
   */
  Scenario();

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
  int coallescing_wells;  //!< the number of identical wells coallescing

  double well_separation;     //!< the separation between coallescing wells
  double dissolution_factor;  //!< the dissolution factor (this should be 1.0!)

  /**
   * @brief Add a stage to the simulation.
   * @param stage the stage definition struct
   * @return the current number of stages
   */
  int add_stage(Stage stage);

  /**
   * @brief Add a debug output to a file
   * @param fout the open steam to print to
   */
  void debug_log(ofstream &fout);
};

/**
 * @brief Class defining the properties for a certain type of rock salt.
 */
class Salt {
 public:
  /**
   * @brief Construct a new Salt object with default vaules for halite,
   * sg_max = 1.2019, rho_s = 2.16 g/cm³. These are default values based on
   * brine created by dissolution of essentially-pure halite.
   */
  Salt(void);

  /**
   * @brief Construct a new Salt object
   * @param sg_max maximum saturated brine density as specific gravity
   * @param rho_solid solid rock density in g/cm³
   */
  Salt(double sg_max, double rho_solid);

  /**
   * @brief Set the saturated sg value
   * @param sg_max maximum saturated brine density, as specific gravity
   */
  void set_sg_max(double sg_max);

  /**
   * @brief Get the saturated sg value
   * @return (double) saturated brine specific gravity
   */
  double get_sg_max(void);

  /**
   * @brief Set the solid density value
   * @param rho_solid solid rock salt density in g/cm³
   */
  void set_solid_density(double rho_solid);

  /**
   * @brief Get the solid density value
   * @return solid rock salt density in g/cm³
   */
  double get_solid_density(void);

  /**
   * @brief Get a COPY of the a-parameter VALUES
   * @return (array<6>) a
   */
  array<double, 6> get_recession_rate_coeff();

  /**
   * @brief Set the recession rate params VALUES
   * @param new_coeff_a length-6 array with parameter values
   */
  void set_recession_rate_coeff(array<double, 6> new_coeff_a);

  /**
   * @brief Get a COPY of the c-parameter VALUES
   * @return array<double,3> c
   */
  array<double, 3> get_sg_wt_pct_convert_coeff();

  /**
   * @brief Set the wt pct sg params VALUES
   * @param new_coeff_c length-3 array
   */
  void set_density_conversion_coeff(array<double, 3> new_coeff_c);

  /**
   * @brief Calculate the weight-percent of brine.
   * @param sg specific gravity of the brine
   * @param temp the brine temperature in degrees Fahrenheit, by default 75
   * @return weight percent of salt in brine
   */
  double get_wt_pct(double sg, double temp = 75.0);

  /**
   * @brief Calculate the specific gravity of brine.
   * @param wt_pct weight-percent salt in the brine
   * @param temp brine temperature in degrees Fahrenheit, by default 75
   * @return specific gravity of the brine
   */
  double get_sg(double wt_pct, double temp = 75.0);

  /**
   * @brief Calculate the recession rate of the wall.
   * @param sg specific gravity of the brine
   * @return recession rate in ft/s
   */
  double get_recession_rate(double sg);

 private:
  double C_sat;  //!< specific gravity of saturated brine
  double rho_s;  //!< density of the surrounding rock in g/cm^3

  array<double, 6> a_;  //!< the coeff for wall recession rate formula
  array<double, 3> c_;  //!< the coeff for converting b/w sg & wt pct
};

/**
 * @brief Model of the water jet flow out the end of tubing.
 */
class JetModel {
 public:
  /**
   * @brief Create a new version-1 jet model.
   */
  JetModel(void);

  /**
   * @brief Create a new jet model, specifying the version.
   * @param ver the version number, 0 for off
   */
  JetModel(int ver);

  /**
   * @brief Set the jet model version
   * @param ver version number, 0 for off
   */
  void set_version(int ver);

  /**
   * @brief Get the jet model version
   * @returns the version
   */
  int get_version(void);

  /**
   * @brief Calculate the jet velocity.
   * @param Q the flow rate in ft3/h
   * @param r the orifice radius in ft
   * @returns the velocity in ft/s
   */
  double velocity(double Q, double r);

  /**
   * @brief Calculate the jet length.
   * @param u the velocity in ft/s
   * @returns the penetration depth in ft
   */
  double length(double u);

 private:
  int version;  //!< the jet model version - 0 is off
};

/**
 * @brief Interface for a class that provides n-derivatives for an ODE func.
 */
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
  /**
   * @brief Create a new, blank ODE solver object.
   */
  ODESolver(void);

  /**
   * @brief Create a new solver object with a specific ODE in mind.
   * @param n the number of equations
   * @param f the derivative function object
   */
  ODESolver(int n, Derivable *f);

  /**
   * @brief Integrate a system of neqn first order ODEs
   * @param y solution vector at t
   * @param t independent variable
   * @param tout point at which solution is desired
   * @param relerr relative error tolerance for local error test
   * @param abserr absolute error tolerance for local error test
   * @param iflag status of integration
   */
  void ode(vector<double> &y, double &t, double tout, double &relerr,
           double &abserr, int &iflag);

  /**
   * @brief Get the number of equations this solver was configured for.
   * @returns number of equations
   */
  int num_eqn(void);

 protected:
  /**
   * @brief Differential equation solver
   * @param y solution vector at T
   * @param t independent variable
   * @param tout point at which solution is desired
   * @param relerr relative error tolerance
   * @param abserr absolute error tolerance
   * @param iflag status of integration
   * @param phase1 work variable
   * @param start work variable
   * @param ns work variable
   * @param nornd work variable
   * @param k work variable
   * @param kold work variable
   * @param isnold work variable
   */
  void de(vector<double> &y, double &t, double tout, double &relerr,
          double &abserr, int &iflag, bool &phase1, bool &start, int &ns,
          bool &nornd, int &k, int &kold, int &isnold);

  /**
   * @brief Approximate the solution near x by polymonial.
   * @param xout points at which solution is known
   * @param yout solution at xout
   * @param kold pass from sansmic::step1 - unchanged
   */
  void intrp(double xout, vector<double> &yout, int kold);

  /**
   * @brief Take a single solver step, used indirectly through ODE
   * @param eps local error tolerance
   * @param start logical variable set true for first step
   * @param k appropriate order for next step
   * @param kold order used for last successful step
   * @param crash logical variable set true when no step can be taken
   * @param phase1 elimainate local retention of variables
   * @param ns elimainate local retention of variables
   * @param nornd elimainate local retention of variables
   */
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

/**
 * @brief Class defining the function for the plume rise integration in ODE.
 */
class PlumeRise : public Derivable {
 public:
  /**
   * @brief Construct a new Plume Rise object
   * @param delta_z cell size
   * @param alpha_coeff entrainment coefficient
   * @param conc concentration vector
   */
  PlumeRise(double delta_z, double alpha_entr, vector<double> &conc);

  /**
   * @brief Calculate derivatives for ODE
   * @param x point to evaluate
   * @param y solution vector
   * @param yp derivative vector
   */
  void func(double &x, vector<double> &y, vector<double> &yp);

  ODESolver solver;    //!< solver for the ode
  const int neqn = 3;  //!< number of equations used

 private:
  double dz;     //!< cell size in ft
  double alpha;  //!< entrainment coefficient

  vector<double> co;  //!< concentration vector for the cavern
};

/**
 * @brief Structure containing the results of the simulation
 */
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

/**
 * @brief The main SANSMIC model class.
 */
class Model {
 public:
  /**
   * @brief Create a new model for the simulation
   * @param out_prefix the prefix to use for file outputs
   */
  Model(string out_prefix = "sansmic");

  int verb = 0;  //!< verbosity setting for output

  /**
   * @brief Configure the model with the scenario provided.
   */
  void configure(Scenario scen);

  /**
   * @brief Run the complete model.
   */
  void run_sim(void);

  /**
   * @brief Run the next stage and all timesteps therein
   * @return 1 if the stage is complete
   */
  int run_stage();

  /**
   * @brief Run the next time step
   * @return 1 if the stage is complete
   */
  int run_step();

  /**
   * @brief Move on to the next stage
   * @return the current stage number
   */
  int init_stage(void);

  /**
   * @brief End the current stage
   * @return 1 if the stage is complete
   */
  int end_stage();

  /**
   * @brief Open output files and initialize the headers.
   * @param append whether the files should be opened with append, by default
   * false
   * @return success code
   */
  int open_outfiles(bool append = false);

  /**
   * @brief Close the file handles for output files
   */
  void close_outfiles(void);

  /**
   * @brief Get the current stage number
   * @return the stage number
   */
  int get_current_stage(void) { return stageNum; }

  /**
   * @brief Get the current time
   * @return the time
   */
  double get_current_time(void) { return days + timet; }

  /**
   * @brief Get the current cavern volume (in bbl)
   * @return the volume
   */
  double get_current_volume(void) { return V_tot; }

  /**
   * @brief Get the single-timestep state of the model in a Results object.
   * @return Results object with single timestep of data.
   */
  Results get_current_state(void);

  /**
   * @brief Get the compelte results object
   * @return Results
   */
  Results get_results(void) { return results; }

  /**
   * @brief Get the stages object
   * @return vector of all stage definitions
   */
  vector<Stage> get_stages(void) { return stages; }

  /**
   * @brief Get the current number of stages
   * @return the number of stages
   */
  size_t num_stages(void) { return stages.size(); }

  /**
   * @brief Is this model currently running?
   * @return running status
   */
  bool get_running_status(void) { return b_running; }

 private:
  string prefix;  //!< output file prefix
  string Q_fill_table_name;
  string Q_inj_table_name;

  ofstream fileOut;  //!< output file
  ofstream fileLog;  //!< log file
  ofstream fileTst;  //!< debug ".tst" file

  Results results;        //!< results object
  Salt salt;              //!< (CSAT,CLAM,DENSAL,CWS) salt properties object
  JetModel jet_model;     //!< the jet model
  PlumeRise *plume_rise;  //!< plume rise model, create on next_stage
  GeomFormat dataFmt;     //!< the format of the geometry data

  vector<Stage> stages;  // simulation stages

  bool b_initialized;     // has the model been initialized
  bool b_running;         // continue running the timestp
  bool b_injecting;       // is the stage injecting or is it static
  bool b_times_up;        // has this step reached the end of the stage
  bool b_obi_below_roof;  // did the obi exit the cavern
  bool b_first_stage;     // is this the first stage or a subsequent stage
  bool b_use_fill_table;  // use a flow table for product fill
  bool b_use_inj_table;   // use a flow table for brine injection
  bool b_just_saved;      // did we just save? don't save twice!

  int n_wells;            // number of coallescing caverns
  int n_cells;            // number of cells
  int n_nodes;            // number of nodes
  int print_freq;         //
  int stepNum;            // current step number
  int stageNum;           // the current stage number
  int obiCell;            // cell containing the interface
  int prodCell;           // production string EOT cell
  int injCell;            // cell containing the injection EOT
  int jetPlumeCell;       // cell containing top of plume
  int noFlowCell;         //
  int startCell;          //
  int injCellBelow;       //
  int obiCellBelow;       //
  int jetPlumeCellBelow;  // index cell below jet
  int maxProdOrJet;       // highest cell that is either production or jet
  int minProdOrJet;       // lowest cell that is either production or jet
  int plumeModelVer;      // the version of the plume model
  int tempModelVer;       // the version of the temperature model
  int im;                 // i minus 1
  int ip;                 // i plus 1
  int L;                  // TODO: ew ew ew change this var name
  int m;                  // mass
  int nEqn = 3;           //
  int idqf;               //
  int idqi;               //
  int i_obi;              // temp blanket cell index
  int i_obiOld;           // temp blanket cell index
  int izbs;               // current OBI cell index
  int jplp;               // bounded jetPlumeCell + 1
  int jpls;               // jet plume cell

  double dt;            //!< timestep size (h)
  double err_cfac;      //!< convergence level
  double C_cavAve;      //!< average cavern brine sg
  double L_jet;         //!< jet model: jet length
  double Q_fOld;        //!< oil fill: last timestep
  double Q_fTot;        //!< oil fill: total volume
  double Q_iOld;        //!< injection: last timestep
  double Q_iTot;        //!< injection: total volume
  double r_inj0;        //!< jet model: injection point radius
  double days;          //!< time: current loop time, in days
  double timet;         //!< time: total time, in days
  double V_ullage;      //!< volume: ullage available
  double u_inj0;        //!< jet model: injection point velocity
  double V_insol;       //!< volume: insolubles in cavern
  double volInsolVent;  //!< volume: insolubles vented out
  double V_usable;      //!< volume: total usable
  double V_tot;         //!< volume: total volume
  double h_obi;         //!< height: OBI/blanket
  double h_insol;       //!< height: of top of insolubles
  double abserr;        //!< ODE solver: absolute tolerance
  double relerr;        //!< ODE solver: relative tolerance
  double floor_depth;   //!< measured depth: bottom of the cavern
  double C_plm;         // plume concentraion
  double ca;            // temp variable A
  double cb;            // temp variable B
  double cc;            // temp variable C
  double cnext;         // next concentration
  double cold;          // old concentration
  double dC_dt;         // concentration derivative w.r.t. time
  double dum;           // temp variable
  double duml;          // temp variable
  double fallf;         // fraction of falling insolubles
  double z_obi;         // obi depth
  double p1, p2, p3, p4, p5, p7, p8, p9, p10;
  double rho_sat;        // saturated density
  double r;              // radius
  double S_d;            //
  double slctim;         //
  double thet;           // wall angle
  double u;              // velocity
  double vdif;           //
  double vel;            //
  double velfac;         //
  double vi;             //
  double vip;            //
  double V_used;         //
  double volRemoved;     //
  double vrn;            //
  double w;              //
  double watsal;         //
  double watsol;         //
  double x;              //
  double dz_inc;         //
  double z_ullage;       //
  double alpha;          // entrainment coefficient
  double V_cell;         //
  double beta;           // diffusion beta coefficient
  double C_bar;          // average specific gravity in plume
  double c1;             //
  double c2;             //
  double c3;             //
  double cpl;            //
  double cpln;           //
  double C_hat;          // specific density saturated brine
  double D_0;            // eddy diffusion coefficient
  double D_mol;          // molecular diffusion coefficient
  double dayOld;         // used for determining daily output
  double C_wall;         // solid density
  double diffCoeff;      //
  double f_dis0;         // base adjustment of dissolution rate from input
  double dt_dz2;         // dt / dz^2
  double dt_min;         // minimum recommended dt based on injection rates
  double dz;             // height of each cell
  double Q_fillBPD;      //
  double Q_fill;         //
  double Q_fillSav;      //
  double Q_inBPD;        //
  double Q_in;           //
  double Q_iSav;         //
  double Q_out;          //
  double Q_outBPD;       //
  double f_insol0;       // fraction insolubles
  double r_outerPipeID;  // inner radius of outer casing
  double r_csgOD;        // outer radius of outer casing
  double dr_dt;          // recession rate
  double r_cavMax;       // maximum cavern diameter
  double r_cavMin;       // minimum cavern diameter
  double r_stag;         // radius of plume at plume stagnation level
  double r_sqd;          // um, volume not radius squared. Actually = r^2 dz
  double r_tbgID;        // inner radius of inner pipe
  double r_innerPipeOD;  // outer radius of inner pipe
  double runMode;        // run mode
  double dx_sep;         // separation between coallescing wells
  double C_cav0;         // initial specific gravity in cavern
  double C_inj;          // specific gravity of injected water
  double sig;            //
  double z_obi_stop;     // OBI-final (for stage termination criteria)
  double ppp;            //
  double stopCriteria;   //
  double t;              // current time
  double t_last;         //
  double t_tot;          //
  double temp;           // temperature
  double t_end;          // time to end injection
  double m_brine;        //
  double m_brineNew;     //
  double m_brineOld;     //
  double m_saltRemove;   //
  double V_plume;        //
  double m_plume;        //
  double t_wait;         // workover period duration
  double vflo;           // flow into plume
  double V_insolRemain;  // insolubles within an cell
  double V_stop;         // stop when oil remaining less than value
  double vppl;           // plume volume
  double vtpl;           // total plume volume height
  double w_d;            // downwind differencing variable
  double w_hat;          // wt pct. of saturated brine
  double w_u;            // upwind differencing variable
  double h_obiSave;      // OBI saved between stages
  double h_inj;          //
  double h_max;          //
  double h_prd;          //
  double h_uso;          //

  vector<int> f_disType;  //!< dissolution regime indicator

  vector<double> akd_prt;       //!< AKD for printing
  vector<double> amd_prt;       //!< AMD for printing
  vector<double> ca_prt;        //!< CA for printing
  vector<double> C_tmp;         //!< concentration
  vector<double> C_cav;         //!< cavern brine sg vector
  vector<double> dC;            //!< change in concentration between steps
  vector<double> z_cav;         //!< cell measured depths
  vector<double> f_dis_prt;     //!< dissolution factor (for output)
  vector<double> dr_prt;        //!< delta radius (for output)
  vector<double> phi;           //!< wall angle of the cell
  vector<double> r_plume;       //!< plume radius vector
  vector<double> C_plume;       //!< plume concentration vector
  vector<double> u_plume;       //!< plume velocity vector
  vector<double> r_cav0;        //!< initial cavern radius
  vector<double> r_cav;         //!< current cavern radius
  vector<double> x_incl;        //!< wall angle factor
  vector<double> h_cav;         //!< cell elevations from floor
  vector<double> tanTheta;      //!< the tangent of the wall angle
  vector<double> V_saltRemove;  //!< volume salt removed
  vector<double> vec_A;         // A-overline in (2.8)
  vector<double> vec_B;         // B-overline in (2.8)
  vector<double> vec_C;         // C-overline in (2.8)
  vector<double> vec_D;         // D-overline in (2.8)
  vector<double> aa;            //
  vector<double> cosTheta;      //
  vector<double> cb_prt;        // used for printing CB
  vector<double> depdkr;        // temp depth vector
  vector<double> f_dis;         // adjusted ratio of dissolution rate
  vector<double> f_disSav;  // base ratio of actual-to-assumed dissolution rate
  vector<double> f_insol;   // volume ratio of insolubles to salt
  vector<double> p2d;
  vector<double> qtim;
  vector<double> qval;
  vector<double> rcscr;
  vector<double> ss;
  vector<double> tunc;
  vector<double> voldkr;       // temp variable for ullage
  vector<double> V_injSigned;  //

  array<int, 2> idxq;   // USED FOR QI/QF TABLES
  array<int, 2> iqnum;  // USED FOR QI/QF TABLES

  /// @brief Perform the leaching
  /// @return is the stage complete
  int leach();

  /// @brief Calculate the brine mass
  /// @param dz the height to calculate within
  /// @return the mass
  double brine_mass(double dz);

  /// @brief Remove salt from the walls
  /// @param i the cell to look at
  /// @param recession_rate the recession rate
  /// @param dt the timestep size
  /// @param dz the vertical cell size
  void remove_salt(int i, double recession_rate, double dt, double dz);

  /// @brief Calculate the plume concentration
  void plume(double ci, double zi, double zb, double &x, double &u, double dz,
             double alpha, int nmax, double &r, double &cpl);

  /// @brief Calculate the new slope
  /// @param i the cell to use
  void slope(int i);

  /// @brief Save results in a results object
  /// @param new_results the object to add to
  /// @param to_file output to files as well
  void save_results(Results &new_results, bool to_file);

  void write_header(ofstream &);
  void write_tst_header(void);
  void write_tst_step(int stage, bool inject);
  void write_tst_end_or_wo(const char *text);
  void write_log_end_stage(void);
  void write_out_timestep_summary(void);
  void write_out_timestep_per_cell(int i, double p1, double p2);
  void write_out_timestep_totals(double qo);
  void write_out_timestep_insolubles(double p1, double p2);
  void write_out_timestep_removed(double p1);
};

/**
 * @brief Namespace for SANSMIC error codes and processing
 */
namespace error {
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

}  // namespace error

}  // namespace sansmic
#endif
