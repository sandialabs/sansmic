// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file util.cpp
 * @brief Utility maths functions
 */

#include <vector>

#include "libsansmic.hpp"

/**
 * @brief Shortcut for square
 * @details Because FORTRAN and Python have an operator for power, a shortcut
 * function was written for brevity and to ensure there are no mistakes due to
 * typos between two complex groupings of products/sums that are then squared.
 * @param val
 * @return double
 */
double sansmic::sq(double val) { return val * val; }

/**
 * @brief absolute value for doubles
 * @param d1 value
 * @return double
 */
double sansmic::dabs(double d1) {
  if (d1 >= 0.0) return d1;
  return -d1;
}

/**
 * @brief implementation of FORTRAN DDIM1 function in C++
 * @param d1 value 1
 * @param d2 value 2
 * @return double
 */
double sansmic::ddim(double d1, double d2) {
  if (d1 > d2) return d1 - d2;
  return 0.0;
}

/**
 * @brief implementation of FORTRAN DSIGN function in C++
 * @param d1
 * @param d2
 * @return double
 */
double sansmic::sign(double d1, double d2) {
  if (d2 >= 0) return sansmic::dabs(d1);
  return -sansmic::dabs(d1);
}

/**
 * @brief implementation of FORTRAN ISIGN function in C++
 * @param i1
 * @param i2
 * @return int
 */
int sansmic::sign(int i1, int i2) {
  if (i2 >= 0) return std::abs(i1);
  return -std::abs(i1);
}

/**
 * @brief interpolation function
 * @param xf x value at which to interpolate
 * @param x known x values
 * @param y known y values
 * @param nmax length of arrays
 * @return double, interpolated y value at xf
 */
double sansmic::xnterp_dkr(double xf, std::vector<double> x,
                           std::vector<double> y, int nmax) {
  double yf = -999.0;
  for (int i = 0; i < nmax; i++) {
    if (xf < x[i] && xf >= x[i + 1]) {
      yf = y[i] + (xf - x[i]) / (x[i + 1] - x[i]) * (y[i + 1] - y[i]);
      break;
    }
  }
  return yf;
}

/**
 * @brief bound a value within a certain range
 * @param value original value
 * @param low bound on low end
 * @param high bound on high end
 * @return double bounded value
 */
double sansmic::bound(double value, double low, double high) {
  return std::min(std::max(value, low), high);
}

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
void sansmic::trigad(int ns, int nf, std::vector<double> &y,
                     std::vector<double> &a, std::vector<double> &b,
                     std::vector<double> &c, std::vector<double> &d) {
  int nsp, nfm, nfd, km, kp;
  double an, bn, dn, anou;
  std::vector<double> cn = std::vector<double>(a.size() + 1, 0.0);

  nsp = ns + 1;
  if (dabs(b[nsp]) <= 1.0e-14) {
    nfm = nf - 1;
    for (int i = nsp; i <= nfm; i++) {
      cn[i + 1] = (d[i] - a[i] * y[i - 1]) / c[i];
      if (i == nsp) {
        y[i] = (y[i + 1] + y[i - 1]) * 0.5;
      }
      y[i + 1] = y[i + 1] - b[i] * y[i] / c[i];
    }
  } else {
    std::vector<double> u(a.size() + 1, 0.0);
    std::vector<double> cntmp(a.size() + 1, 0.0);
    u[ns] = 1.0;
    cn[ns] = 0.0;
    nfd = nf - nsp;
    for (int k = nsp; k <= nf; k++) {
      km = k - 1;
      an = a[k];
      bn = b[k];
      cn[k] = c[k];
      dn = d[k];
      anou = an / u[km];
      u[k] = bn - anou * cn[km];
      y[k] = dn - anou * y[km];
    }
    y[nf] = y[nf] - y[nf + 1] * cn[nf];
    y[nf] = y[nf] / u[nf];
    for (int kk = 1; kk <= nfd; kk++) {
      int k = nf - kk;
      kp = k + 1;
      y[k] = (y[k] - cn[k] * y[kp]) / u[k];
    }
  }

  return;
}
