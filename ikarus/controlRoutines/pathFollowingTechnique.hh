// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <memory>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include "ikarus/solver/nonLinearSolver/newtonRaphson.hh"
#include <ikarus/controlRoutines/adaptiveStepSizing.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphsonWithScalarSubsidiaryFunction.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observerMessages.hh>
#include <ikarus/utils/pathFollowingFunctions.hh>

namespace Ikarus {

  struct PathFollowingInformation {
    bool success{false};
  };

  template <typename NonLinearOperator>
  double getMinEigenValue(NonLinearOperator& nonOp, const int eigenValueOfInterest, bool updateDisplacements = false) {
    const auto& K = nonOp.derivative();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;  // Dense matrices only
    es.compute(K);
    auto eigenValues  = es.eigenvalues();
    auto eigenVectors = es.eigenvectors();
    if (updateDisplacements) nonOp.firstParameter() = eigenVectors.col(eigenValueOfInterest);
    return eigenValues[eigenValueOfInterest];  // only works because eigenValues here is sorted already by Eigen
  }

  template <typename ScalarType = double>
  struct BisectionMethod {
    template <typename NonLinearOperator, typename NonLinearSolver, typename PathFollowingType, typename PathFollowing>
    void operator()(NonLinearOperator& nonOp, const std::shared_ptr<NonLinearSolver>& nonLinearSolver,
                    SubsidiaryArgs& subsidiaryArgs, PathFollowingType& pathFollowingType, SolverInfos& solverInfo,
                    const int eigenValueOfInterest, PathFollowing& pf) {
      std::cout << "########################################" << std::endl;
      std::cout << "######  BISECTION METHOD BEGINS  #######" << std::endl;
      std::cout << "########################################" << std::endl;
      auto minEigenValue         = getMinEigenValue(nonOp, eigenValueOfInterest);
      int numberOfBisections     = 100;
      int bisectionCounter       = 0;
      auto subsidiaryFunctionPtr = [&](auto&& args) { return pathFollowingType.evaluateSubsidiaryFunction(args); };

      /// Bisection Method
      while ((abs(minEigenValue) > 1e-8) and (bisectionCounter < numberOfBisections)) {
        pathFollowingType.intermediatePredictionBisectionForward(nonOp, subsidiaryArgs);
        solverInfo    = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);
        minEigenValue = getMinEigenValue(nonOp, eigenValueOfInterest);
        bisectionCounter++;
        if (minEigenValue < 0.0) pathFollowingType.intermediatePredictionBack(nonOp, subsidiaryArgs);
        std::cout << std::setprecision(16) << "minEigenValue:\t" << minEigenValue << "\tbisectionCounter:\t"
                  << bisectionCounter << std::endl;
      }
      pf.notify(ControlMessages::CRITICAL_POINT_CHANGED);
      minEigenValue = getMinEigenValue(nonOp, eigenValueOfInterest, true);
    }
  };

  template <typename NonLinearSolver, typename PathFollowingType = Ikarus::StandardArcLength,
            typename AdaptiveStepSizing     = Ikarus::DefaultAdaptiveStepSizing<>,
            typename CriticalPointEstimator = Ikarus::BisectionMethod<>>
  requires Concepts::PathFollowingStrategy<PathFollowingType, typename NonLinearSolver::NonLinearOperator>
  class PathFollowing : public IObservable<ControlMessages> {
  public:
    PathFollowing(const std::shared_ptr<NonLinearSolver>& p_nonLinearSolver, int loadSteps, double stepSize,
                  PathFollowingType p_pathFollowingType   = Ikarus::StandardArcLength{},
                  AdaptiveStepSizing p_adaptiveStepSizing = {}, bool p_determineCriticalPoint = false,
                  CriticalPointEstimator p_criticalPointEstimator = {})
        : nonLinearSolver{p_nonLinearSolver},
          loadSteps_{loadSteps},
          stepSize_{stepSize},
          pathFollowingType_{p_pathFollowingType},
          adaptiveStepSizing{p_adaptiveStepSizing},
          determineCriticalPoint{p_determineCriticalPoint},
          criticalPointEstimator{p_criticalPointEstimator} {}

    PathFollowingInformation run() {
      spdlog::info("Started path following with subsidiary equation: {}", pathFollowingType_.name);
      PathFollowingInformation info({false});
      auto& nonOp = nonLinearSolver->nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED);

      SubsidiaryArgs subsidiaryArgs;

      subsidiaryArgs.stepSize  = stepSize_;
      subsidiaryArgs.loadSteps = loadSteps_;
      subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
      subsidiaryArgs.DD.setZero();
      /// adaptive step-sizing
      if (std::is_same_v<Ikarus::FancyAdaptiveStepSizing<>, decltype(adaptiveStepSizing)>) {
        auto numEigen = 3;
        subsidiaryArgs.eigenValues.setZero(numEigen, subsidiaryArgs.loadSteps);
      }

      /// Initializing solver
      this->notify(ControlMessages::STEP_STARTED);
      auto subsidiaryFunctionPtr = [&](auto&& args) { return pathFollowingType_.evaluateSubsidiaryFunction(args); };
      pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs);
      auto solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);
      if (not solverInfo.success) return info;
      this->notify(ControlMessages::SOLUTION_CHANGED);
      this->notify(ControlMessages::STEP_ENDED);
      int bisectionCounter       = 0;
      int eigenValueOfInterest   = 0;
      int numberOfCriticalPoints = 4;

      /// Calculate predictor for a particular step
      for (int ls = 0; ls < loadSteps_; ++ls) {
        this->notify(ControlMessages::STEP_STARTED);
        subsidiaryArgs.actStep = ls;
        if (ls > 0) pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);

        solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);

        auto& DDummy1      = nonOp.firstParameter();
        auto& LambdaDummy1 = nonOp.lastParameter();

        auto DDummy2      = DDummy1;
        auto LambdaDummy2 = LambdaDummy1;

        if (determineCriticalPoint and (getMinEigenValue(nonOp, eigenValueOfInterest) < 0.0)) {
          pathFollowingType_.intermediatePredictionBack(nonOp, subsidiaryArgs);
          criticalPointEstimator(nonOp, nonLinearSolver, subsidiaryArgs, pathFollowingType_, solverInfo,
                                 eigenValueOfInterest, *this);
          this->notify(ControlMessages::CRITICAL_POINT_CHANGED);
          bisectionCounter++;
          eigenValueOfInterest++;
          std::cout << "########################################" << std::endl;
          std::cout << "######   BISECTION METHOD ENDS   #######" << std::endl;
          std::cout << "########################################" << std::endl;
          if (bisectionCounter == numberOfCriticalPoints) {
            std::cout << "####### bisectionCounter reached the maximum value #######" << std::endl;
            break;
          }
        }
        DDummy1      = DDummy2;
        LambdaDummy1 = LambdaDummy2;

        adaptiveStepSizing(solverInfo, nonOp, subsidiaryArgs);

        //
        if (not solverInfo.success) return info;
        this->notify(ControlMessages::SOLUTION_CHANGED);
        this->notify(ControlMessages::STEP_ENDED, subsidiaryArgs.stepSize);
      }
      this->notify(ControlMessages::CONTROL_ENDED);
      info.success = true;
      return info;
    }

  private:
    std::shared_ptr<NonLinearSolver> nonLinearSolver;
    int loadSteps_;
    double stepSize_;
    PathFollowingType pathFollowingType_;
    AdaptiveStepSizing adaptiveStepSizing;
    CriticalPointEstimator criticalPointEstimator;
    bool determineCriticalPoint;
  };
}  // namespace Ikarus
