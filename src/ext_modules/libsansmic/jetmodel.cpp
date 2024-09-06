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
 * @brief Create a new version-1 jet model.
 */
sansmic::JetModel::JetModel(void) { version = 1; }

/**
 * @brief Create a new jet model, specifying the version.
 * @param ver the version number, 0 for off
 */
sansmic::JetModel::JetModel(int ver) { version = ver; }

/**
 * @brief Set the jet model version
 * @param ver version number, 0 for off
 */
void sansmic::JetModel::set_version(int ver) { version = ver; }

/**
 * @brief Get the jet model version
 * @returns the version
 */
int sansmic::JetModel::get_version(void) { return version; }

/**
 * @brief Calculate the jet velocity.
 * @param Q the flow rate in ft3/h
 * @param r the orifice radius in ft
 * @returns the velocity in ft/s
 */
double sansmic::JetModel::velocity(double Q, double r) {
  /*                    Q [ft3/h]
   *  u_jet0 = -----------------------------
   *            3600 [s/h] _pi_ r^2 [ft2]
   */
  return Q / (hour * _pi_ * sq(r));
}

/**
 * @brief Calculate the jet length.
 * @param u the velocity in ft/s
 * @returns the penetration depth in ft
 */
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
