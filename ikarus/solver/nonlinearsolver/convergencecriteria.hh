// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file convergencecriteria.hh
 * \brief A collection of convergence criteria used by the nonlinear solvers.
 */

#include <dune/common/float_cmp.hh>

#include <ikarus/utils/linearalgebrahelper.hh>

#pragma once

namespace Ikarus::ConvergenceCriteria {

struct ResiduumNorm
{
  ResiduumNorm() = default;
  /**
   * \brief Call operator when the norm of the residuum is chosen as the convergence criteria.
   * \tparam NLO Type of the nonlinear operator.
   * \tparam SS Type of the nonlinear solver settings.
   * \tparam CV Type of the correction vector used to update the solution.
   * \param nonLinearOperator The nonlinear operator.
   * \param solverSettings The nonlinear solver settings.
   * \param corectionVector The correction vector used to update the solution.
   */
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const std::shared_ptr<NLO>& nonLinearOperator, const NSS& solverSettings,
                  [[maybe_unused]] const CV& correctionVector) const {
    static_assert(std::is_same_v<CV, std::remove_cvref_t<decltype(nonLinearOperator->firstParameter())>>,
                  "CV type does not match nonLinearOperator->firstParameter()");
    const auto& rx = nonLinearOperator->value();
    auto rNorm     = norm(rx);
    return Dune::FloatCmp::le(static_cast<double>(rNorm), solverSettings.tol);
  }
};

struct CorrectionNorm
{
  CorrectionNorm() = default;
  /**
   * \brief Call operator when the norm of the correction vector is chosen as the convergence criteria.
   * \tparam NLO Type of the nonlinear operator.
   * \tparam SS Type of the nonlinear solver settings.
   * \tparam CV Type of the correction vector used to update the solution.
   * \param nonLinearOperator The nonlinear operator.
   * \param solverSettings The nonlinear solver settings.
   * \param corectionVector The correction vector used to update the solution.
   */
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const std::shared_ptr<NLO>& nonLinearOperator, const NSS& solverSettings,
                  [[maybe_unused]] const CV& correctionVector) const {
    static_assert(std::is_same_v<CV, std::remove_cvref_t<decltype(nonLinearOperator->firstParameter())>>,
                  "CV type does not match nonLinearOperator->firstParameter()");
    auto dNorm = norm(correctionVector);
    return Dune::FloatCmp::le(static_cast<double>(dNorm), solverSettings.tol);
  }
};

struct MaximumSolution
{
  /**
   * \brief Constructor for MaximumSolution.
   * \param maxSol The expected maximum value of the solution vector.
   */
  explicit MaximumSolution(double maxSol)
      : maxSol_{maxSol} {};
  /**
   * \brief Call operator when the maximum value of the solution vector is chosen as the convergence criteria.
   * \tparam NLO Type of the nonlinear operator.
   * \tparam SS Type of the nonlinear solver settings.
   * \tparam CV Type of the correction vector used to update the solution.
   * \param nonLinearOperator The nonlinear operator.
   * \param solverSettings The nonlinear solver settings.
   * \param corectionVector The correction vector used to update the solution.
   */
  template <typename NLO, typename NSS, typename CV>
  bool operator()(const std::shared_ptr<NLO>& nonLinearOperator, const NSS& solverSettings,
                  [[maybe_unused]] const CV& correctionVector) const {
    static_assert(std::is_same_v<CV, std::remove_cvref_t<decltype(nonLinearOperator->firstParameter())>>,
                  "CV type does not match nonLinearOperator->firstParameter()");
    const auto& sol = nonLinearOperator->firstParameter();
    auto maxSol     = max(sol);
    return Dune::FloatCmp::le(static_cast<double>(std::abs(maxSol - maxSol_)), solverSettings.tol);
  }

private:
  double maxSol_;
};

} // namespace Ikarus::ConvergenceCriteria
