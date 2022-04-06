//
// Created by lex on 30/12/2021.
//

#pragma once
#include <iostream>

#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <ikarus/utils/observer/observerMessages.hh>
#include <ikarus/linearAlgebra/linearAlgebraHelper.hh>
#include <ikarus/utils/observer/observer.hh>

namespace Ikarus {

  struct NonlinearSolverSettings {
    double tol{1e-8};
    int maxIter{20};
  };

  struct SolverInformation {
    bool sucess{false};
    double residualnorm{0.0};
    int iterations{0};
  };

  template <typename LinearSolver, typename MatrixType, typename VectorType>
  concept LinearSolverC = requires(LinearSolver& linearSolver, MatrixType& Ax, VectorType& vec) {
    linearSolver.analyzePattern(Ax);
    linearSolver.factorize(Ax);
    linearSolver.solve(vec);
  };

  template <typename NonLinearOperatorImpl,
            typename LinearSolver = std::function<typename NonLinearOperatorImpl::ValueType(
                const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>,
            typename UpdateType
            = std::conditional_t<std::is_floating_point_v<typename NonLinearOperatorImpl::template Parameter<0>>,
                                 typename NonLinearOperatorImpl::template Parameter<0>, Eigen::VectorXd>>
  class NewtonRaphson : public IObservable<NonLinearSolverMessages> {
  public:
    using LinearSolverScalarFunctionType = std::function<typename NonLinearOperatorImpl::ValueType(
        const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>;

    static constexpr bool isLinearSolver = LinearSolverC<LinearSolver, typename NonLinearOperatorImpl::DerivativeType,
                                                         typename NonLinearOperatorImpl::ValueType>;

    using ResultType         = typename NonLinearOperatorImpl::template Parameter<0>;
    using UpdateFunctionType = std::function<void(ResultType&, const UpdateType&)>;

    explicit NewtonRaphson(
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
          updateFunction{p_updateFunction} {
      if constexpr (std::is_same_v<typename NonLinearOperatorImpl::ValueType, Eigen::VectorXd>)
        corr.setZero(nonLinearOperator().value().size());
    }

    using NonLinearOperator = NonLinearOperatorImpl;

    void setup(const NonlinearSolverSettings& p_settings) { settings = p_settings; }

    struct NoPredictor {};
    template <typename SolutionType = NoPredictor>
    requires std::is_same_v<SolutionType, NoPredictor> || std::is_convertible_v<
        SolutionType, std::remove_cvref_t<typename NonLinearOperatorImpl::ValueType>>
        SolverInformation solve(const SolutionType& dx_predictor = NoPredictor{}) {
      this->notify(NonLinearSolverMessages::INIT);
      SolverInformation solverInformation;
      solverInformation.sucess = true;
      nonLinearOperator().updateAll();
      const auto& rx = nonLinearOperator().value();
      const auto& Ax = nonLinearOperator().derivative();
      auto& x        = nonLinearOperator().firstParameter();
      auto rNorm     = norm(rx);
      decltype(rNorm) dNorm;
      int iter{0};
      if constexpr (not std::is_same_v<SolutionType, NoPredictor>) updateFunction(x, dx_predictor);
      if constexpr (isLinearSolver) linearSolver.analyzePattern(Ax);
      while (rNorm > settings.tol && iter <= settings.maxIter) {
        this->notify(NonLinearSolverMessages::ITERATION_STARTED);
        if constexpr (isLinearSolver) {
          linearSolver.factorize(Ax);
          corr  = -linearSolver.solve(rx);
          dNorm = corr.norm();
          updateFunction(x, corr);
        } else {
          corr  = -linearSolver(rx, Ax);
          dNorm = norm(corr);
          updateFunction(x, corr);
        }
        this->notify(NonLinearSolverMessages::CORRECTIONNORM_UPDATED, dNorm);
        this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
        nonLinearOperator().updateAll();
        rNorm = norm(rx);
        this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, rNorm);
        this->notify(NonLinearSolverMessages::ITERATION_ENDED);
        ++iter;
      }
      solverInformation.iterations   = iter;
      solverInformation.residualnorm = rNorm;
      this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, iter, rNorm, settings.tol);
      return solverInformation;
    }

    auto& nonLinearOperator() { return nonLinearOperator_; }

  private:
    NonLinearOperatorImpl nonLinearOperator_;
    typename NonLinearOperatorImpl::ValueType corr;
    LinearSolver linearSolver;
    UpdateFunctionType updateFunction;
    NonlinearSolverSettings settings;
  };

  template <typename NonLinearOperatorImpl,
            typename LinearSolver = std::function<typename NonLinearOperatorImpl::ValueType(
                const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>,
            typename UpdateType
            = std::conditional_t<std::is_floating_point_v<typename NonLinearOperatorImpl::template Parameter<0>>,
                                 typename NonLinearOperatorImpl::template Parameter<0>, Eigen::VectorXd>>
 auto makeNewtonRaphson(
      const NonLinearOperatorImpl& p_nonLinearOperator,
      LinearSolver&& p_linearSolver = [](const typename NonLinearOperatorImpl::ValueType& a,
                                         const typename NonLinearOperatorImpl::ValueType& b) { return a / b; },
      std::function<void(typename NonLinearOperatorImpl::template Parameter<0>&, const UpdateType&)> p_updateFunction =
          [](typename NonLinearOperatorImpl::template Parameter<0>& a, const UpdateType& b) {
            using Ikarus::operator+=;
            a += b;
          }) {
    return std::make_shared<NewtonRaphson<NonLinearOperatorImpl, LinearSolver, UpdateType>>(
        p_nonLinearOperator, std::move(p_linearSolver), std::move(p_updateFunction));
  }

}  // namespace Ikarus