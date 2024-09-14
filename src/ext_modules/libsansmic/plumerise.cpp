// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file plumerise.cpp
 * @brief Plume rise derivative class
 */

#include <cmath>
#include <vector>

#include "libsansmic.hpp"

using namespace sansmic;

/**
 * @brief Construct a new Plume Rise object
 * @param delta_z cell size
 * @param alpha_coeff entrainment coefficient
 * @param conc concentration vector
 */
PlumeRise::PlumeRise(double delta_z, double alpha_entr,
                     std::vector<double> &conc) {
  dz = delta_z;
  alpha = alpha_entr;
  co = conc;
  solver = ODESolver(neqn, (Derivable *)this);
}

/**
 * @brief Calculate derivatives for ODE
 * @param x point to evaluate
 * @param y solution vector
 * @param yp derivative vector
 */
void PlumeRise::func(double &x, std::vector<double> &y,
                     std::vector<double> &yp) {
  double u = 1.0e-32;
  double cp, gg;
  if (std::abs(y[1]) > 1.0e-6) {
    u = y[2] / y[1];
  }
  int j = int(x / dz) + 1;
  cp = (co[j + 1] - co[j]) / dz;
  gg = std_gravity / foot * cp;
  yp[1] = 2.0 * alpha * std::sqrt(std::abs(y[2]));
  yp[2] = 2.0 * y[3] / u;
  yp[3] = 2.0 * y[1] * gg;
}
