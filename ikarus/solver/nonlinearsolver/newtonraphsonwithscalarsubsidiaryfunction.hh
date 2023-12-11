// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <iosfwd>
#include <utility>

#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>
#include <ikarus/utils/pathfollowingfunctions.hh>

namespace Ikarus {

  struct NewtonRaphsonWithSubsidiaryFunctionSettings {
    double tol{1e-8};
    int maxIter{30};
  };

  struct SolverInfos {
    bool success{false};
    double residualnorm{0.0};
    int iterations{0};
  };

  template <typename LinearSolver, typename MatrixType, typename VectorType>
  concept LinearSolverCheck = requires(LinearSolver& linearSolver, MatrixType& Ax, VectorType& vec) {
    linearSolver.analyzePattern(Ax);
    linearSolver.factorize(Ax);
    linearSolver.solve(vec, vec);
  };

  template <typename NonLinearOperatorImpl,
            typename LinearSolver = std::function<typename NonLinearOperatorImpl::ValueType(
                const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>,
            typename UpdateType
            = std::conditional_t<std::is_floating_point_v<typename NonLinearOperatorImpl::template ParameterValue<0>>,
                                 typename NonLinearOperatorImpl::template ParameterValue<0>, Eigen::VectorXd>>
  class NewtonRaphsonWithSubsidiaryFunction : public IObservable<NonLinearSolverMessages> {
  public:
    using LinearSolverScalarFunctionType = std::function<typename NonLinearOperatorImpl::ValueType(
        const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>;

    static constexpr bool isLinearSolver
        = LinearSolverCheck<LinearSolver, typename NonLinearOperatorImpl::DerivativeType,
                            typename NonLinearOperatorImpl::ValueType>;

    using ResultType         = typename NonLinearOperatorImpl::template ParameterValue<0>;
    using UpdateFunctionType = std::function<void(ResultType&, const UpdateType&)>;

    explicit NewtonRaphsonWithSubsidiaryFunction(
        const NonLinearOperatorImpl& p_nonLinearOperator,
        LinearSolver&& p_linearSolver = [](const typename NonLinearOperatorImpl::ValueType& a,
                                           const typename NonLinearOperatorImpl::ValueType& b) { return a / b; },
        std::function<void(ResultType&, const UpdateType&)> p_updateFunction =
            [](ResultType& a, const UpdateType& b) {
              using Ikarus::operator+=;
              a += b;
            })
        : nonLinearOperator_{p_nonLinearOperator},
          linearSolver{std::move(p_linearSolver)},
          updateFunction{p_updateFunction} {}

    using NonLinearOperator = NonLinearOperatorImpl;

    void setup(const NewtonRaphsonWithSubsidiaryFunctionSettings& p_settings) { settings = p_settings; }

    struct NoPredictor {};
    template <typename SolutionType = NoPredictor, typename SubsidiaryType>
    requires std::is_same_v<SolutionType, NoPredictor> || std::is_convertible_v<
        SolutionType, std::remove_cvref_t<typename NonLinearOperatorImpl::ValueType>>
        SolverInfos solve(SubsidiaryType& subsidiaryFunction, SubsidiaryArgs& subsidiaryArgs) {
      this->notify(NonLinearSolverMessages::INIT);

      /// Initializations
      SolverInfos solverInformation;
      solverInformation.success = true;
      auto& x                   = nonLinearOperator().firstParameter();  // x = D (Displacements)
      auto& lambda              = nonLinearOperator().lastParameter();

      /// Determine Fext0
      /// It is assumed that Fext = Fext0 * lambda such that dRdlambda = Fext0
      /// Generalization for Fext0 = Fext0(lambda) is not implemented

      auto lambdaDummy = lambda;
      auto DDummy      = x;
      x.setZero();
      lambda = 1.0;
      nonLinearOperator().template update<0>();
      const auto Fext0 = (-nonLinearOperator().value()).eval();  // dRdlambda(lambda)
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
      if constexpr (isLinearSolver) linearSolver.analyzePattern(Ax);

      /// Iterative solving scheme
      while (rNorm > settings.tol && iter < settings.maxIter) {
        this->notify(NonLinearSolverMessages::ITERATION_STARTED);

        /// Two-step solving procedure
        residual2d.resize(rx.rows(), 2);
        residual2d << -rx, Fext0;
        sol2d.resize(rx.rows(), 2);

        if constexpr (isLinearSolver) {
          linearSolver.factorize(Ax);
          linearSolver.solve(sol2d, residual2d);
        } else {
          sol2d = linearSolver(residual2d, Ax);
        }

        subsidiaryFunction(subsidiaryArgs);
        this->notify(NonLinearSolverMessages::SCALARSUBSIDIARY_UPDATED, subsidiaryArgs.f);

        const double deltalambda = (-subsidiaryArgs.f - subsidiaryArgs.dfdDD.dot(sol2d.col(0)))
                                   / (subsidiaryArgs.dfdDD.dot(sol2d.col(1)) + subsidiaryArgs.dfdDlambda);
        deltaD = sol2d.col(0) + deltalambda * sol2d.col(1);

        updateFunction(x, deltaD);
        updateFunction(subsidiaryArgs.DD, deltaD);

        lambda += deltalambda;
        subsidiaryArgs.Dlambda += deltalambda;

        dNorm = sqrt(deltaD.dot(deltaD) + deltalambda * deltalambda);
        nonLinearOperator().updateAll();
        rNorm = sqrt(rx.dot(rx) + subsidiaryArgs.f * subsidiaryArgs.f);

        this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
        this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, dNorm);
        this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, rNorm);
        this->notify(NonLinearSolverMessages::ITERATION_ENDED);
        ++iter;
      }

      if (iter == settings.maxIter) solverInformation.success = false;
      solverInformation.iterations   = iter;
      solverInformation.residualnorm = rNorm;
      if (solverInformation.success)
        this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, iter, rNorm, settings.tol);

      return solverInformation;
    }

    auto& nonLinearOperator() { return nonLinearOperator_; }

  private:
    NonLinearOperatorImpl nonLinearOperator_;
    LinearSolver linearSolver;
    UpdateFunctionType updateFunction;
    NewtonRaphsonWithSubsidiaryFunctionSettings settings;
  };

  template <typename NonLinearOperatorImpl,
            typename LinearSolver = std::function<typename NonLinearOperatorImpl::ValueType(
                const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>,
            typename UpdateType
            = std::conditional_t<std::is_floating_point_v<typename NonLinearOperatorImpl::template ParameterValue<0>>,
                                 typename NonLinearOperatorImpl::template ParameterValue<0>, Eigen::VectorXd>>
  auto makeNewtonRaphsonWithSubsidiaryFunction(
      const NonLinearOperatorImpl& p_nonLinearOperator,
      LinearSolver&& p_linearSolver = [](const typename NonLinearOperatorImpl::ValueType& a,
                                         const typename NonLinearOperatorImpl::ValueType& b) { return a / b; },
      std::function<void(typename NonLinearOperatorImpl::template ParameterValue<0>&, const UpdateType&)>
          p_updateFunction
      =
          [](typename NonLinearOperatorImpl::template ParameterValue<0>& a, const UpdateType& b) {
            using Ikarus::operator+=;
            a += b;
          }) {
    return std::make_shared<NewtonRaphsonWithSubsidiaryFunction<NonLinearOperatorImpl, LinearSolver, UpdateType>>(
        p_nonLinearOperator, std::forward<LinearSolver>(p_linearSolver), std::move(p_updateFunction));
  }
}  // namespace Ikarus
