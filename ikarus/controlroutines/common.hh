// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file common.hh
 * \brief Common things from the control routines.
 */

#pragma once

#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

/**
 * \brief A helper function that returns a linear solver suitable for symmetric, positive-definite sparse or dense
 * matrices, based on the nonlinear solver.
 *
 * \tparam NLS Type of the nonlinear solver.
 *
 * \param nls The nonlinear solver.
 *
 * \return A LinearSolver.
 */
template <typename NLS>
auto createSPDLinearSolverFromNonLinearSolver(const NLS& nls) {
  using JacobianType = std::remove_cvref_t<typename NLS::JacobianType>;
  static_assert((traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, JacobianType>::value) or
                    (traits::isSpecializationTypeNonTypeAndType<Eigen::SparseMatrix, JacobianType>::value),
                "Linear solver not implemented for the chosen derivative type of the non-linear operator");
  SolverTypeTag solverTag;

  auto&& residual = nls.residual();

  if constexpr (traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, typename NLS::JacobianType>::value)
    solverTag = SolverTypeTag::d_LDLT;
  else
    solverTag = SolverTypeTag::sd_SimplicialLDLT;
  return LinearSolver(solverTag);
}

template <typename NLS>
typename NLS::Domain::SolutionVectorType idbcIncrement(typename NLS::Domain& x, const NLS& nls, double Dlambda) {
  auto y = x;
  y.parameter() += Dlambda;
  y.syncParameterAndGlobalSolution(nls.updateFunction());
  const auto delta = (y.globalSolution() - x.globalSolution()).eval();
  return delta;
}
} // namespace Ikarus
