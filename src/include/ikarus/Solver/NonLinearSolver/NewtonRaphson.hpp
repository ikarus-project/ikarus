//
// Created by lex on 30/12/2021.
//

#pragma once
#include <iostream>

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/Solver/LinearSolver/LinearSolver.h"
#include "ikarus/utils/LinearAlgebraHelper.h"
#include "ikarus/utils/Observer/observerMessages.h"
#include <ikarus/utils/Observer/observer.h>

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
  template <typename NonLinearOperatorImpl,
            typename LinearSolver   = std::function<typename NonLinearOperatorImpl::ValueType(
                  const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>,
            typename UpdateFunction = std::function<void(typename NonLinearOperatorImpl::ValueType&,
                                                         const typename NonLinearOperatorImpl::ValueType&)> >
  class NewtonRaphson : public IObservable<NonLinearSolverMessages> {
  public:
    using LinearSolverStdFunctionType = std::function<typename NonLinearOperatorImpl::ValueType(
        const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>;

    explicit NewtonRaphson(
        const NonLinearOperatorImpl& p_nonLinearOperator,
        LinearSolver&& p_linearSolver   = [](const typename NonLinearOperatorImpl::ValueType& a,
                                           const typename NonLinearOperatorImpl::ValueType& b) { return a / b; },
        UpdateFunction p_updateFunction = [](typename NonLinearOperatorImpl::ValueType& a,
                                             const typename NonLinearOperatorImpl::ValueType& b) { a += b; })
        : nonLinearOperator_{p_nonLinearOperator},
          linearSolver{std::move(p_linearSolver)},
          updateFunction{p_updateFunction} {}

    using NonLinearOperator = NonLinearOperatorImpl;

    void setup(const NonlinearSolverSettings& p_settings) { settings = p_settings; }

    struct NoPredictor {};
    template <typename SolutionType = NoPredictor>
      requires std::is_same_v<SolutionType,
                              NoPredictor> || std::is_same_v<SolutionType, typename NonLinearOperatorImpl::ValueType>
    SolverInformation solve(const SolutionType& dx_predictor = NoPredictor{}) {
      this->notify(NonLinearSolverMessages::INIT);
      SolverInformation solverInformation;
      solverInformation.sucess = true;
      nonLinearOperator().updateAll();
      const auto& rx = nonLinearOperator().value();
      const auto& Ax = nonLinearOperator().derivative();
      auto& x        = nonLinearOperator().template nthParameter<0>();
      auto rNorm     = norm(rx);
      decltype(rNorm) dNorm;
      int iter{0};
      if constexpr (not std::is_same_v<SolutionType, NoPredictor>) updateFunction(x, dx_predictor);
      if constexpr (!std::is_same_v<LinearSolver, LinearSolverStdFunctionType>) linearSolver.analyzePattern(Ax);
      while (rNorm > settings.tol && iter <= settings.maxIter) {
        this->notify(NonLinearSolverMessages::ITERATION_STARTED);
        if constexpr (!std::is_same_v<LinearSolver, LinearSolverStdFunctionType>) {
          linearSolver.factorize(Ax);
          const Eigen::VectorXd D = -linearSolver.solve(rx);
          dNorm                   = D.norm();
          updateFunction(x, D);
        } else {
          const auto D = -linearSolver(rx, Ax);
          dNorm        = D;
          updateFunction(x, D);
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
    LinearSolver linearSolver;
    UpdateFunction updateFunction;
    NonlinearSolverSettings settings;
  };

}  // namespace Ikarus