// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file jetmodel.cpp
 * @brief Jet model implementation
 */

#include "libsansmic.hpp"

using namespace sansmic;

/**
 * @brief Create the default jet model, which is version 1.
 */
sansmic::JetModel::JetModel(void) { version = 1; }

/**
 * @brief Create a jet model with a specific version.
 * @param ver the version to use
 */
sansmic::JetModel::JetModel(int ver) { version = ver; }

void sansmic::JetModel::set_version(int ver) { version = ver; }

int sansmic::JetModel::get_version(void) { return version; }

double sansmic::JetModel::velocity(double Q, double r) {
  /*                    Q [ft3/h]
   *  u_jet0 = -----------------------------
   *            3600 [s/h] _pi_ r^2 [ft2]
   */
  return Q / (hour * _pi_ * sq(r));
}

double sansmic::JetModel::length(double u) {
  if (version == 0) {
    return 0.0;
  }
  /*
   *  l_jet = 0.5 [s] * u_inj0 [ft/s]
   *
   */
  return 0.5 * u;
}
