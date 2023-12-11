// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <iosfwd>

#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/utils/defaultfunctions.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

  struct NewtonRaphsonSettings {
    double tol{1e-8};
    int maxIter{20};
  };

  struct SolverInformation {
    explicit operator bool() const { return success; }
    bool success{false};
    double residualnorm{0.0};
    int iterations{0};
  };

  template <typename LinearSolver, typename MatrixType, typename VectorType>
  concept LinearSolverC = requires(LinearSolver& linearSolver, MatrixType& Ax, VectorType& vec) {
    linearSolver.analyzePattern(Ax);
    linearSolver.factorize(Ax);
    linearSolver.solve(vec, vec);
  };

  template <typename NonLinearOperatorImpl, typename LinearSolver = SolverDefault,
            typename UpdateFunctionType_ = UpdateDefault>
  class NewtonRaphson : public IObservable<NonLinearSolverMessages> {
  public:
    static constexpr bool isLinearSolver = LinearSolverC<LinearSolver, typename NonLinearOperatorImpl::DerivativeType,
                                                         typename NonLinearOperatorImpl::ValueType>;

    using ResultType         = typename NonLinearOperatorImpl::template ParameterValue<0>;
    using UpdateFunctionType = UpdateFunctionType_;

    explicit NewtonRaphson(const NonLinearOperatorImpl& p_nonLinearOperator, LinearSolver&& p_linearSolver = {},
                           UpdateFunctionType_ p_updateFunction = {})
        : nonLinearOperator_{p_nonLinearOperator},
          linearSolver{std::move(p_linearSolver)},
          updateFunction{p_updateFunction} {
      if constexpr (std::is_same_v<typename NonLinearOperatorImpl::ValueType, Eigen::VectorXd>)
        correction.setZero(nonLinearOperator().value().size());
    }

    using NonLinearOperator = NonLinearOperatorImpl;

    void setup(const NewtonRaphsonSettings& p_settings) { settings = p_settings; }

    struct NoPredictor {};
    template <typename SolutionType = NoPredictor>
    requires std::is_same_v<SolutionType, NoPredictor> || std::is_convertible_v<
        SolutionType, std::remove_cvref_t<typename NonLinearOperatorImpl::ValueType>>
    [[nodiscard(
        "The solve method returns information of the solution process. You should store this information and check if "
        "it was successful")]] SolverInformation
    solve(const SolutionType& dx_predictor = NoPredictor{}) {
      this->notify(NonLinearSolverMessages::INIT);
      SolverInformation solverInformation;
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
      solverInformation.iterations   = iter;
      solverInformation.residualnorm = static_cast<double>(rNorm);
      if (solverInformation.success)
        this->notify(NonLinearSolverMessages::FINISHED_SUCESSFULLY, iter, static_cast<double>(rNorm), settings.tol);
      return solverInformation;
    }

    auto& nonLinearOperator() { return nonLinearOperator_; }

  private:
    NonLinearOperatorImpl nonLinearOperator_;
    typename NonLinearOperatorImpl::ValueType correction;
    LinearSolver linearSolver;
    UpdateFunctionType updateFunction;
    NewtonRaphsonSettings settings;
  };

  template <typename NonLinearOperatorImpl, typename LinearSolver = SolverDefault, typename Update = UpdateDefault>
  auto makeNewtonRaphson(const NonLinearOperatorImpl& p_nonLinearOperator, LinearSolver&& p_linearSolver = {},
                         Update&& p_updateFunction = {}) {
    return std::make_shared<NewtonRaphson<NonLinearOperatorImpl, LinearSolver, Update>>(
        p_nonLinearOperator, std::forward<LinearSolver>(p_linearSolver), std::move(p_updateFunction));
  }

}  // namespace Ikarus
