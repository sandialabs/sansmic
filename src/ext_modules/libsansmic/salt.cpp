// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file salt.cpp
 * @brief Salt class implementation
 */

#include <array>
#include <cmath>

#include "libsansmic.hpp"

/**
 * @brief Construct a new Salt object with default vaules for halite.
 */
sansmic::Salt::Salt(void) {
  C_sat = 1.2019;  // sg
  rho_s = 2.16;    // g/cm^3
  a_[0] = 7.6091661297e-1;
  a_[1] = -3.8715516531;
  a_[2] = 7.8254117638;
  a_[3] = -7.8395924067;
  a_[4] = 3.8789476581;
  a_[5] = -0.7533873462;
  c_[0] = 1.0019795678;
  c_[1] = 0.67459890538;
  c_[2] = 0.32350531044;
}

/**
 * @brief Construct a new Salt object
 * @param sg_max maximum saturated brine density as specific gravity
 * @param rho_solid solid rock density in g/cm³
 */
sansmic::Salt::Salt(double sg_max, double rho_solid) {
  C_sat = sg_max;
  rho_s = rho_solid;
  a_[0] = 7.6091661297e-1;
  a_[1] = -3.8715516531;
  a_[2] = 7.8254117638;
  a_[3] = -7.8395924067;
  a_[4] = 3.8789476581;
  a_[5] = -0.7533873462;
  c_[0] = 1.0019795678;
  c_[1] = 0.67459890538;
  c_[2] = 0.32350531044;
}

/**
 * @brief Set the saturated sg value
 * @param sg_max maximum saturated brine density, as specific gravity
 */
void sansmic::Salt::set_sg_max(double sg_max) { C_sat = sg_max; }

/**
 * @brief Get the saturated sg value
 * @return (double) saturated brine specific gravity
 */
double sansmic::Salt::get_sg_max(void) { return C_sat; }

/**
 * @brief Set the solid density value
 * @param rho_solid solid rock salt density in g/cm³
 */
void sansmic::Salt::set_solid_density(double rho_solid) { rho_s = rho_solid; }

/**
 * @brief Get the solid density value
 * @return solid rock salt density in g/cm³
 */
double sansmic::Salt::get_solid_density(void) { return rho_s; }

/**
 * @brief Set the recession rate params VALUES
 * @param new_coeff_a length-6 array with parameter values
 */
void sansmic::Salt::set_recession_rate_coeff(
    std::array<double, 6> new_coeff_a) {
  for (int i = 0; i < 6; i++) {
    a_[i] = new_coeff_a[i];
  }
}

/**
 * @brief Set the wt pct sg params VALUES
 * @param new_coeff_c length-3 array
 */
void sansmic::Salt::set_density_conversion_coeff(
    std::array<double, 3> new_coeff_c) {
  for (int i = 0; i < 3; i++) {
    c_[i] = new_coeff_c[i];
  }
}

/**
 * @brief Get a COPY of the a-parameter VALUES
 * @return (array<6>) a
 */
std::array<double, 6> sansmic::Salt::get_recession_rate_coeff() { return a_; }

/**
 * @brief Get a COPY of the c-parameter VALUES
 * @return array<double,3> c
 */
std::array<double, 3> sansmic::Salt::get_sg_wt_pct_convert_coeff() {
  return c_;
}

/**
 * @brief Calculate the weight-percent of brine.
 * @param sg specific gravity of the brine
 * @param temp the brine temperature in degrees Fahrenheit, by default 75
 * @return weight percent of salt in brine
 */
double sansmic::Salt::wt_pct(double sg, double temp) {
  double x, wt_pct, delta;
  if (sg <= 1.0726) {
    x = (sg - 1.0) / 0.726;
  } else if (sg >= C_sat) {
    x = 0.2632;
  } else {
    delta = c_[1] * c_[1] + 4.0 * c_[2] * (sg - c_[0]);
    x = 0.5 * (sqrt(delta) - c_[1]) / c_[2];
  }
  wt_pct = x * (1.0 + 0.000353 * (temp - 75.0));
  return wt_pct;
}

/**
 * @brief Calculate the specific gravity of brine.
 * @param wt_pct weight-percent salt in the brine
 * @param temp brine temperature in degrees Fahrenheit, by default 75
 * @return specific gravity of the brine
 */
double sansmic::Salt::sg(double wt_pct, double temp) {
  double y, sg;
  y = wt_pct / (1.0 + 0.000353 * (temp - 75.0));
  if (y <= 0.1) {
    sg = 1.0 + 0.726 * y;
  } else if (y >= 0.2632) {
    sg = C_sat;
  } else {
    sg = fmin(c_[0] + c_[1] * y + c_[2] * y * y, C_sat);
  }
  return sg;
}

/**
 * @brief Calculate the recession rate of the wall.
 * @param sg specific gravity of the brine
 * @return recession rate in ft/s
 */
double sansmic::Salt::recession_rate(double sg) {
  double z = 0;
  if (sg == 0.0) {
    z = 1.0;
  } else if (sg >= C_sat) {
    z = 0.0;
  } else {
    z = a_[0];
    for (int k = 1; k < 6; k++) {
      z = sg * z + a_[k];
    }
    z = z / sg;
    if (sg >= 1.2019) z = 0;
  };
  return z;
}
