// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file defaultfunctions.hh
 * \brief Collection of fallback default functions
 */

#pragma once
#include "linearalgebrahelper.hh"

#include <type_traits>

namespace Ikarus::utils {
/**
 * \struct SolverDefault
 * \ingroup utils
 * \brief Default functor for solving operations.
 * \details This functor provides a default implementation for solving operations by performing division.
 * It is intended to be used in generic contexts where a default solver is needed.
 * \tparam A Type of the left operand.
 * \tparam B Type of the right operand.
 * \param a The left operand.
 * \param b The right operand.
 * \return The result of the division operation (a / b).
 */
struct SolverDefault
{
  template <typename A, typename B>
  constexpr auto operator()(A&& a, B&& b) const {
    return a / b;
  }
};

/**
 * \struct SyncFERequirements
 * \brief A helper struct that tells the update function of the nonlinear solver to sync the global solution based on
 * the parameter, i.e., to sync the FERequirements.
 */
struct SyncFERequirements
{
} inline constexpr syncFERequirements;

/**
 * \struct IDBCForceDefault
 * \brief Default struct used to represent that no inhomogeneous Dirichlet BCs are present.
 */
struct IDBCForceDefault
{
};

/**
 * \struct UpdateDefault
 * \ingroup utils
 * \brief Default functor for updating operations.
 * \details This functor provides a default implementation for updating operations by performing addition.
 * It is intended to be used in generic contexts where a default updater is needed.
 * \tparam A Type of the target operand to be updated.
 * \tparam B Type of the value to be added for updating.
 * \param a The target operand to be updated.
 * \param b The value to be added for updating.
 * \return The result of the addition operation (a += b).
 */
struct UpdateDefault
{
  template <typename A, typename B = SyncFERequirements>
  constexpr void operator()(A&& a, B&& b = {}) const {
    if constexpr (not std::is_same_v<B, SyncFERequirements>) {
      using Ikarus::operator+=;
      a += b;
    }
  }
};

} // namespace Ikarus::utils
