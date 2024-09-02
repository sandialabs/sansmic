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

double sansmic::sq(double val) { return val * val; }

double sansmic::dabs(double d1) {
  if (d1 >= 0.0) return d1;
  return -d1;
}

double sansmic::ddim(double d1, double d2) {
  if (d1 > d2) return d1 - d2;
  return 0.0;
}

double sansmic::sign(double d1, double d2) {
  if (d2 >= 0) return sansmic::dabs(d1);
  return -sansmic::dabs(d1);
}

int sansmic::sign(int i1, int i2) {
  if (i2 >= 0) return std::abs(i1);
  return -std::abs(i1);
}

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

double sansmic::bound(double value, double low, double high) {
  return std::min(std::max(value, low), high);
}

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
