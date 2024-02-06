// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file newtonraphson.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

/**
 * \struct NewtonRaphsonSettings
 * \brief Settings for the Newton-Raphson solver.
 */
struct NewtonRaphsonSettings
{
  double tol{1e-8};
  int maxIter{20};
};

/**
 * \class NewtonRaphson
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \relates makeNewtonRaphson
 * \ingroup solvers
 */
template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
class NewtonRaphson : public IObservable<NonLinearSolverMessages>
{
public:
  ///< Compile-time boolean indicating if the linear solver satisfies the non-linear solver concept
  static constexpr bool isLinearSolver =
      Ikarus::Concepts::LinearSolverCheck<LinearSolver, typename NLO::DerivativeType, typename NLO::ValueType>;

  ///< Type representing the parameter vector of the nonlinear operator.
  using ValueType = typename NLO::template ParameterValue<0>;

  using UpdateFunction    = UF;  ///< Type representing the update function.
  using NonLinearOperator = NLO; ///< Type of the non-linear operator

  /**
   * \brief Constructor for NewtonRaphson.
   * \param nonLinearOperator Nonlinear operator to solve.
   * \param linearSolver Linear solver used internally (default is SolverDefault).
   * \param updateFunction Update function (default is UpdateDefault).
   */
  explicit NewtonRaphson(const NonLinearOperator& nonLinearOperator, LS&& linearSolver = {},
                         UpdateFunction updateFunction = {})
      : nonLinearOperator_{nonLinearOperator},
        linearSolver_{std::move(linearSolver)},
        updateFunction_{updateFunction} {
    if constexpr (std::is_same_v<typename NonLinearOperator::ValueType, Eigen::VectorXd>)
      correction_.setZero(this->nonLinearOperator().value().size());
  }

  /**
   * \brief Set up the solver with the given settings.
   * \param settings Newton-Raphson settings.
   */
  void setup(const NewtonRaphsonSettings& settings) { settings_ = settings; }

#ifndef DOXYGEN
  struct NoPredictor
  {
  };
#endif
  /**
   * \brief Solve the nonlinear system.
   * \param dxPredictor Predictor for the solution increment (default is NoPredictor).
   * \return Information about the solution process.
   */
  template <typename SolutionType = NoPredictor>
  requires std::is_same_v<SolutionType, NoPredictor> ||
           std::is_convertible_v<SolutionType, std::remove_cvref_t<typename NonLinearOperator::ValueType>>
  [[nodiscard(
      "The solve method returns information of the solution process. You should store this information and check if "
      "it was successful")]] Ikarus::NonLinearSolverInformation
  solve(const SolutionType& dxPredictor = NoPredictor{}) {
    this->notify(NonLinearSolverMessages::INIT);
    Ikarus::NonLinearSolverInformation solverInformation;
    solverInformation.success = true;
    auto& x                   = nonLinearOperator().firstParameter();
    if constexpr (not std::is_same_v<SolutionType, NoPredictor>)
      updateFunction_(x, dxPredictor);
    nonLinearOperator().updateAll();
    const auto& rx = nonLinearOperator().value();
    const auto& Ax = nonLinearOperator().derivative();
    auto rNorm     = norm(rx);
    decltype(rNorm) dNorm;
    int iter{0};
    if constexpr (isLinearSolver)
      linearSolver_.analyzePattern(Ax);
    while (rNorm > settings_.tol && iter < settings_.maxIter) {
      this->notify(NonLinearSolverMessages::ITERATION_STARTED);
      if constexpr (isLinearSolver) {
        linearSolver_.factorize(Ax);
        linearSolver_.solve(correction_, -rx);
        dNorm = correction_.norm();
        updateFunction_(x, correction_);
      } else {
        correction_ = -linearSolver_(rx, Ax);
        dNorm       = norm(correction_);
        updateFunction_(x, correction_);
      }
      this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, static_cast<double>(dNorm));
      this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
      nonLinearOperator().updateAll();
      rNorm = norm(rx);
      this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, static_cast<double>(rNorm));
      this->notify(NonLinearSolverMessages::ITERATION_ENDED);
      ++iter;
    }
    if (iter == settings_.maxIter)
      solverInformation.success = false;
    solverInformation.iterations     = iter;
    solverInformation.residualNorm   = static_cast<double>(rNorm);
    solverInformation.correctionNorm = static_cast<double>(dNorm);
    if (solverInformation.success)
      this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, iter);
    return solverInformation;
  }

  /**
   * \brief Access the nonlinear operator.
   * \return Reference to the nonlinear operator.
   */
  auto& nonLinearOperator() { return nonLinearOperator_; }

private:
  NonLinearOperator nonLinearOperator_;
  typename NonLinearOperator::ValueType correction_;
  LS linearSolver_;
  UpdateFunction updateFunction_;
  NewtonRaphsonSettings settings_;
};

/**
 * \brief Function to create a NewtonRaphson solver instance.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param nonLinearOperator Nonlinear operator to solve.
 * \param linearSolver Linear solver used internally (default is SolverDefault).
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
auto makeNewtonRaphson(const NLO& nonLinearOperator, LS&& linearSolver = {}, UF&& updateFunction = {}) {
  return std::make_shared<NewtonRaphson<NLO, LS, UF>>(nonLinearOperator, std::forward<LS>(linearSolver),
                                                      std::move(updateFunction));
}

} // namespace Ikarus
