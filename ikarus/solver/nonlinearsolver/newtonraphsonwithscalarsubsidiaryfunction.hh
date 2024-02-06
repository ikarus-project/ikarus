// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file newtonraphson.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

#include <iosfwd>
#include <utility>

#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {
/**
 * \struct NewtonRaphsonWithSubsidiaryFunctionSettings
 * \brief Settings for the Newton-Raphson solver with subsidiary function.
 */
struct NewtonRaphsonWithSubsidiaryFunctionSettings
{
  double tol{1e-8};
  int maxIter{20};
};

/**
 * \brief Newton-Raphson solver with subsidiary function.
 *
 * This class provides a Newton-Raphson solver for solving nonlinear systems with a subsidiary function.
 * It uses a linear solver to handle the linear system arising in each iteration.
 *
 * \tparam NLO Type of the nonlinear operator.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 */
template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
class NewtonRaphsonWithSubsidiaryFunction : public IObservable<NonLinearSolverMessages>
{
public:
  ///< Compile-time boolean indicating if the linear solver satisfies the non-linear solver concept
  static constexpr bool isLinearSolver =
      Ikarus::Concepts::LinearSolverCheck<LS, typename NLO::DerivativeType, typename NLO::ValueType>;

  ///< Type representing the parameter vector of the nonlinear operator.
  using ValueType = typename NLO::template ParameterValue<0>;
  ///< Type representing the update function.
  using UpdateFunctionType = UF;
  using NonLinearOperator  = NLO; ///< Type of the non-linear operator

  /**
   * \brief Constructor for NewtonRaphsonWithSubsidiaryFunction.
   *
   * \param nonLinearOperator Nonlinear operator to solve.
   * \param linearSolver Linear solver used internally (default is SolverDefault).
   * \param updateFunction Update function (default is UpdateDefault).
   */
  explicit NewtonRaphsonWithSubsidiaryFunction(const NLO& nonLinearOperator, LS&& linearSolver = {},
                                               UpdateFunctionType updateFunction = {})
      : nonLinearOperator_{nonLinearOperator},
        linearSolver_{std::move(linearSolver)},
        updateFunction_{updateFunction} {}

  /**
   * \brief Setup the Newton-Raphson solver with subsidiary function.
   *
   * \param p_settings Settings for the solver.
   */
  void setup(const NewtonRaphsonWithSubsidiaryFunctionSettings& settings) { settings_ = settings; }

#ifndef DOXYGEN
  struct NoPredictor
  {
  };
#endif

  /**
   * \brief Solve the nonlinear system using the Newton-Raphson method with subsidiary function.
   *
   * \tparam SolutionType Type of the solution predictor (default is NoPredictor).
   * \tparam SubsidiaryType Type of the subsidiary function.
   * \param subsidiaryFunction Subsidiary function to be solved.
   * \param subsidiaryArgs Additional arguments for the subsidiary function.
   * \param dxPredictor Predictor for the solution increment (default is NoPredictor).
   * \return Information about the solution process.
   */
  template <typename SolutionType = NoPredictor, typename SubsidiaryType>
  requires std::is_same_v<SolutionType, NoPredictor> ||
           std::is_convertible_v<SolutionType, std::remove_cvref_t<typename NLO::ValueType>>
  [[nodiscard(
      "The solve method returns information of the solution process. You should store this information and check if "
      "it was successful")]] NonLinearSolverInformation
  solve(SubsidiaryType& subsidiaryFunction, SubsidiaryArgs& subsidiaryArgs,
        const SolutionType& dxPredictor = NoPredictor{}) {
    this->notify(NonLinearSolverMessages::INIT);

    /// Initializations
    Ikarus::NonLinearSolverInformation solverInformation;
    solverInformation.success = true;
    auto& x                   = nonLinearOperator().firstParameter(); // x = D (Displacements)
    if constexpr (not std::is_same_v<SolutionType, NoPredictor>)
      updateFunction_(x, dxPredictor);
    auto& lambda = nonLinearOperator().lastParameter();

    /// Determine Fext0
    /// It is assumed that Fext = Fext0 * lambda such that dRdlambda = Fext0
    /// Generalization for Fext0 = Fext0(lambda) is not implemented

    auto lambdaDummy = lambda;
    auto DDummy      = x;
    x.setZero();
    lambda = 1.0;
    nonLinearOperator().template update<0>();
    const auto Fext0 = (-nonLinearOperator().value()).eval(); // dRdlambda(lambda)
    lambda           = lambdaDummy;
    x                = DDummy;

    Eigen::MatrixX2<double> residual2d, sol2d;
    nonLinearOperator().updateAll();
    const auto& rx = nonLinearOperator().value();
    const auto& Ax = nonLinearOperator().derivative();

    Eigen::VectorXd deltaD;
    deltaD.resizeLike(rx);
    deltaD.setZero();

    subsidiaryArgs.dfdDD.resizeLike(Fext0);

    subsidiaryFunction(subsidiaryArgs);
    auto rNorm = sqrt(rx.dot(rx));
    decltype(rNorm) dNorm;
    int iter{0};
    if constexpr (isLinearSolver)
      linearSolver_.analyzePattern(Ax);

    /// Iterative solving scheme
    while (rNorm > settings_.tol && iter < settings_.maxIter) {
      this->notify(NonLinearSolverMessages::ITERATION_STARTED);

      /// Two-step solving procedure
      residual2d.resize(rx.rows(), 2);
      residual2d << -rx, Fext0;
      sol2d.resize(rx.rows(), 2);

      if constexpr (isLinearSolver) {
        linearSolver_.factorize(Ax);
        linearSolver_.solve(sol2d, residual2d);
      } else {
        sol2d = linearSolver(residual2d, Ax);
      }

      subsidiaryFunction(subsidiaryArgs);

      const double deltalambda = (-subsidiaryArgs.f - subsidiaryArgs.dfdDD.dot(sol2d.col(0))) /
                                 (subsidiaryArgs.dfdDD.dot(sol2d.col(1)) + subsidiaryArgs.dfdDlambda);
      deltaD = sol2d.col(0) + deltalambda * sol2d.col(1);

      updateFunction_(x, deltaD);
      updateFunction_(subsidiaryArgs.DD, deltaD);

      lambda += deltalambda;
      subsidiaryArgs.Dlambda += deltalambda;

      dNorm = sqrt(deltaD.dot(deltaD) + deltalambda * deltalambda);
      nonLinearOperator().updateAll();
      rNorm = sqrt(rx.dot(rx) + subsidiaryArgs.f * subsidiaryArgs.f);

      this->notify(NonLinearSolverMessages::SOLUTION_CHANGED, static_cast<double>(lambda));
      this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, static_cast<double>(dNorm));
      this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, static_cast<double>(rNorm));
      this->notify(NonLinearSolverMessages::ITERATION_ENDED);

      ++iter;
    }

    if (iter == settings_.maxIter)
      solverInformation.success = false;
    solverInformation.iterations     = iter;
    solverInformation.residualNorm   = rNorm;
    solverInformation.correctionNorm = dNorm;
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
  NLO nonLinearOperator_;
  LinearSolver linearSolver_;
  UpdateFunctionType updateFunction_;
  NewtonRaphsonWithSubsidiaryFunctionSettings settings_;
};

/**
 * \brief Function to create a NewtonRaphson with subsidiary function solver instance.
 * \tparam NLO Type of the nonlinear operator to solve.
 * \tparam LS Type of the linear solver used internally (default is SolverDefault).
 * \tparam UF Type of the update function (default is UpdateDefault).
 * \param nonLinearOperator Nonlinear operator to solve.
 * \param linearSolver Linear solver used internally (default is SolverDefault).
 * \param updateFunction Update function (default is UpdateDefault).
 * \return Shared pointer to the NewtonRaphson solver instance.
 */
template <typename NLO, typename LS = utils::SolverDefault, typename UF = utils::UpdateDefault>
auto makeNewtonRaphsonWithSubsidiaryFunction(const NLO& nonLinearOperator, LS&& linearSolver = {},
                                             UF&& updateFunction = {}) {
  return std::make_shared<NewtonRaphsonWithSubsidiaryFunction<NLO, LS, UF>>(
      nonLinearOperator, std::forward<LS>(linearSolver), std::move(updateFunction));
}
} // namespace Ikarus
