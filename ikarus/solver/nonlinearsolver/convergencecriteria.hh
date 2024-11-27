// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file convergencecriteria.hh
 * \brief A collection of convergence criteria used by the nonlinear solvers.
 */

#include <dune/common/float_cmp.hh>

#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>

#pragma once

namespace Ikarus::ConvergenceCriteria {

struct ResiduumNorm
{
  /**
   * \brief Call operator when the norm of the residuum is chosen as the convergence criteria.
   * \tparam NLO Type of the nonlinear operator.
   * \tparam SS Type of the nonlinear solver settings.
   * \tparam CV Type of the correction vector used to update the solution.
   * \param nonLinearOperator The nonlinear operator.
   * \param solverSettings The nonlinear solver settings.
   * \param corectionVector The correction vector used to update the solution.
   * \param subsidiaryArgs The subsidiary arguments used with \ref NewtonRaphsonWithSubsidiaryFunction (optional with
   * its value (f=0))
   */
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const NLO& nonLinearOperator, const NSS& solverSettings, const CV& correctionVector,
                  const SubsidiaryArgs& subsidiaryArgs = {.f = 0.0}) {
    const auto& rx = nonLinearOperator.value();
    auto rNorm     = sqrt(rx.dot(rx) + subsidiaryArgs.f * subsidiaryArgs.f);
    return Dune::FloatCmp::le(static_cast<double>(rNorm), solverSettings.tol);
  }
};

struct CorrectionNorm
{
  /**
   * \brief Call operator when the norm of the correction vector is chosen as the convergence criteria.
   * \tparam NLO Type of the nonlinear operator.
   * \tparam SS Type of the nonlinear solver settings.
   * \tparam CV Type of the correction vector used to update the solution.
   * \param nonLinearOperator The nonlinear operator.
   * \param solverSettings The nonlinear solver settings.
   * \param corectionVector The correction vector used to update the solution.
   * \param subsidiaryArgs The subsidiary arguments used with \ref NewtonRaphsonWithSubsidiaryFunction (optional with
   * its value (f=0))
   */
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const NLO& nonLinearOperator, const NSS& solverSettings, const CV& correctionVector,
                  const SubsidiaryArgs& subsidiaryArgs = {.f = 0.0}) {
    auto dNorm = norm(correctionVector);
    return Dune::FloatCmp::le(static_cast<double>(dNorm), solverSettings.tol);
  }
};

struct MaximumCorrection
{
  /**
   * \brief Call operator when the maximum value of the correction vector is chosen as the convergence criteria.
   * \tparam NLO Type of the nonlinear operator.
   * \tparam SS Type of the nonlinear solver settings.
   * \tparam CV Type of the correction vector used to update the solution.
   * \param nonLinearOperator The nonlinear operator.
   * \param solverSettings The nonlinear solver settings.
   * \param corectionVector The correction vector used to update the solution.
   * \param subsidiaryArgs The subsidiary arguments used with \ref NewtonRaphsonWithSubsidiaryFunction (optional with
   * its value (f=0))
   */
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const NLO& nonLinearOperator, const NSS& solverSettings, const CV& correctionVector,
                  const SubsidiaryArgs& subsidiaryArgs = {.f = 0.0}) {
    auto maxCorr = max(correctionVector);
    return Dune::FloatCmp::le(static_cast<double>(maxCorr), solverSettings.tol);
  }
};

} // namespace Ikarus::ConvergenceCriteria
