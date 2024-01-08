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
  struct NewtonRaphsonSettings {
    double tol{1e-8};
    int maxIter{20};
  };

  /**
   * \class NewtonRaphson
   * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
   * \tparam NonLinearOperatorImpl Type of the nonlinear operator to solve.
   * \tparam LinearSolver Type of the linear solver used internally (default is SolverDefault).
   * \tparam UpdateFunctionTypeImpl Type of the update function (default is UpdateDefault).
   * \relates makeNewtonRaphson
   * \ingroup solvers
   */
  template <typename NonLinearOperatorImpl, typename LinearSolver = utils::SolverDefault,
            typename UpdateFunctionTypeImpl = utils::UpdateDefault>
  class NewtonRaphson : public IObservable<NonLinearSolverMessages> {
  public:
    ///< Compile-time boolean indicating if the linear solver satisfies the non-linear solver concept
    static constexpr bool isLinearSolver
        = Ikarus::Concepts::LinearSolverCheck<LinearSolver, typename NonLinearOperatorImpl::DerivativeType,
                                              typename NonLinearOperatorImpl::ValueType>;

    ///< Type representing the parameter vector of the nonlinear operator.
    using ValueType = typename NonLinearOperatorImpl::template ParameterValue<0>;

    using UpdateFunctionType = UpdateFunctionTypeImpl;  ///< Type representing the update function.
    using NonLinearOperator  = NonLinearOperatorImpl;   ///< Type of the non-linear operator

    /**
     * \brief Constructor for NewtonRaphson.
     * \param p_nonLinearOperator Nonlinear operator to solve.
     * \param p_linearSolver Linear solver used internally (default is SolverDefault).
     * \param p_updateFunction Update function (default is UpdateDefault).
     */
    explicit NewtonRaphson(const NonLinearOperatorImpl& p_nonLinearOperator, LinearSolver&& p_linearSolver = {},
                           UpdateFunctionTypeImpl p_updateFunction = {})
        : nonLinearOperator_{p_nonLinearOperator},
          linearSolver{std::move(p_linearSolver)},
          updateFunction{p_updateFunction} {
      if constexpr (std::is_same_v<typename NonLinearOperatorImpl::ValueType, Eigen::VectorXd>)
        correction.setZero(nonLinearOperator().value().size());
    }

    /**
     * \brief Set up the solver with the given settings.
     * \param p_settings Newton-Raphson settings.
     */
    void setup(const NewtonRaphsonSettings& p_settings) { settings = p_settings; }

#ifndef DOXYGEN
    struct NoPredictor {};
#endif
    /**
     * \brief Solve the nonlinear system.
     * \param dx_predictor Predictor for the solution increment (default is NoPredictor).
     * \return Information about the solution process.
     */
    template <typename SolutionType = NoPredictor>
    requires std::is_same_v<SolutionType, NoPredictor> || std::is_convertible_v<
        SolutionType, std::remove_cvref_t<typename NonLinearOperatorImpl::ValueType>>
    [[nodiscard(
        "The solve method returns information of the solution process. You should store this information and check if "
        "it was successful")]] Ikarus::NonLinearSolverInformation
    solve(const SolutionType& dx_predictor = NoPredictor{}) {
      this->notify(NonLinearSolverMessages::INIT);
      Ikarus::NonLinearSolverInformation solverInformation;
      solverInformation.success = true;
      auto& x                   = nonLinearOperator().firstParameter();
      if constexpr (not std::is_same_v<SolutionType, NoPredictor>) updateFunction(x, dx_predictor);
      nonLinearOperator().updateAll();
      const auto& rx = nonLinearOperator().value();
      const auto& Ax = nonLinearOperator().derivative();
      auto rNorm     = norm(rx);
      decltype(rNorm) dNorm;
      int iter{0};
      if constexpr (isLinearSolver) linearSolver.analyzePattern(Ax);
      while (rNorm > settings.tol && iter < settings.maxIter) {
        this->notify(NonLinearSolverMessages::ITERATION_STARTED);
        if constexpr (isLinearSolver) {
          linearSolver.factorize(Ax);
          linearSolver.solve(correction, -rx);
          dNorm = correction.norm();
          updateFunction(x, correction);
        } else {
          correction = -linearSolver(rx, Ax);
          dNorm      = norm(correction);
          updateFunction(x, correction);
        }
        this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, static_cast<double>(dNorm));
        this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
        nonLinearOperator().updateAll();
        rNorm = norm(rx);
        this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, static_cast<double>(rNorm));
        this->notify(NonLinearSolverMessages::ITERATION_ENDED);
        ++iter;
      }
      if (iter == settings.maxIter) solverInformation.success = false;
      solverInformation.iterations     = iter;
      solverInformation.residualNorm   = static_cast<double>(rNorm);
      solverInformation.correctionNorm = static_cast<double>(dNorm);
      if (solverInformation.success) this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, iter);
      return solverInformation;
    }

    /**
     * \brief Access the nonlinear operator.
     * \return Reference to the nonlinear operator.
     */
    auto& nonLinearOperator() { return nonLinearOperator_; }

  private:
    NonLinearOperatorImpl nonLinearOperator_;
    typename NonLinearOperatorImpl::ValueType correction;
    LinearSolver linearSolver;
    UpdateFunctionType updateFunction;
    NewtonRaphsonSettings settings;
  };

  /**
   * \brief Function to create a NewtonRaphson solver instance.
   * \tparam NonLinearOperatorImpl Type of the nonlinear operator to solve.
   * \tparam LinearSolver Type of the linear solver used internally (default is SolverDefault).
   * \tparam UpdateFunctionType Type of the update function (default is UpdateDefault).
   * \param p_nonLinearOperator Nonlinear operator to solve.
   * \param p_linearSolver Linear solver used internally (default is SolverDefault).
   * \param p_updateFunction Update function (default is UpdateDefault).
   * \return Shared pointer to the NewtonRaphson solver instance.
   */
  template <typename NonLinearOperatorImpl, typename LinearSolver = utils::SolverDefault,
            typename UpdateFunctionType = utils::UpdateDefault>
  auto makeNewtonRaphson(const NonLinearOperatorImpl& p_nonLinearOperator, LinearSolver&& p_linearSolver = {},
                         UpdateFunctionType&& p_updateFunction = {}) {
    return std::make_shared<NewtonRaphson<NonLinearOperatorImpl, LinearSolver, UpdateFunctionType>>(
        p_nonLinearOperator, std::forward<LinearSolver>(p_linearSolver), std::move(p_updateFunction));
  }

}  // namespace Ikarus
