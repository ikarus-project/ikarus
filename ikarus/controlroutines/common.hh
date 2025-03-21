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
[[nodiscard]] Eigen::VectorXd predictorForNewLoadLevel(const NLS& nls, const typename NLS::Domain& x_old,
                                                       typename NLS::Domain& x_new) {
  auto&& residual = nls.residual();

  auto&& R = residual(x_new);
  auto&& K = derivative(residual)(x_new);

  x_old.syncParameterAndGlobalSolution(nls.updateFunction());
  x_new.syncParameterAndGlobalSolution(nls.updateFunction());
  auto deltaDFromIHBC = x_new.globalSolution() - x_old.globalSolution();
  R -= K * deltaDFromIHBC; // compute the internal forces due to the displcament increment
                           // F_Int_dir = K* delta_u_dir, delta_u_dir is only non-zero for the inhomogeneous part

  auto linearSolver = createSPDLinearSolverFromNonLinearSolver(nls);
  linearSolver.analyzePattern(K);
  linearSolver.factorize(K);
  Eigen::VectorXd dPredictor;
  linearSolver.solve(dPredictor, -R);

  return dPredictor;
}
} // namespace Ikarus
