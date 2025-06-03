// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file common.hh
 * \brief Common things from the control routines.
 */

#pragma once

#include <concepts>

#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

namespace Concepts {
  /**
   * \concept HasValidIDBCForceFunction
   * \brief A concept to check if the underlying solver has a valid function to handle inhomogeneous Dirichlet BCs.
   * \tparam NLS Type of the nonlinear solver.
   */
  template <typename NLS>
  concept HasValidIDBCForceFunction = not std::same_as<typename NLS::IDBCForceFunction, utils::IDBCForceDefault>;
} // namespace Concepts

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

  if constexpr (traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, JacobianType>::value)
    solverTag = SolverTypeTag::d_LDLT;
  else
    solverTag = SolverTypeTag::sd_SimplicialLDLT;
  return LinearSolver(solverTag);
}

/**
 * \brief A helper function to calculate the increment in the solution vector based on inhomogeneous Dirichlet BCs.
 * \tparam NLS Type of the nonlinear solver.
 * \param x The solution.
 * \param nls The nonlinear solver.
 * \param Dlambda The step size denoting increment in the load factor.
 */
template <typename NLS>
typename NLS::Domain::SolutionVectorType idbcIncrement(typename NLS::Domain& x, const NLS& nls, double Dlambda) {
  if constexpr (Concepts::HasValidIDBCForceFunction<NLS>) {
    auto y = x;
    y.parameter() += Dlambda;
    y.syncParameterAndGlobalSolution(nls.updateFunction());
    const auto delta = (y.globalSolution() - x.globalSolution()).eval();
    return delta;
  } else {
    Eigen::VectorXd v;
    v.setZero(x.globalSolution().size());
    return v;
  }
}
} // namespace Ikarus
