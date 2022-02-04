//
// Created by lex on 17/12/2021.
//

#pragma once
#include <memory>

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/Variables/ParameterFactory.h"
#include <ikarus/utils/Observer/observer.h>
#include <ikarus/utils/Observer/observerMessages.h>

namespace Ikarus {

  template <typename NonLinearSolver,
//            typename FEManager, typename DirichletManager,
            typename... LinearAlgebraFunctionArgs
            >
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(
//        FEManager& feManager,
//                DirichletManager& dirichletManager,
                NonLinearSolver&& p_nonLinearSolver,
                const int loadSteps, const std::pair<double, double>& tbeginEnd)
        :
//          feManager_{&feManager},
//          dirichletManager_{&dirichletManager},
          nonLinearSolver{std::move(p_nonLinearSolver)},
          loadSteps_{loadSteps},
          parameterBegin_{tbeginEnd.first},
          parameterEnd_{tbeginEnd.second},
          stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {}

    void run() {
      auto& nonOp = nonLinearSolver.nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED);
      auto& loadParameter    = nonOp.template nthParameter<0>();
      auto& x    = nonOp.template nthParameter<1>();
      //assemble u
//      auto& PredictRHS    = nonOp.value()(u=(0,0,0,0,uhat));

      loadParameter = 0.0;
//      auto& x                = feManager_->getVariables();
//      x+=u=(0,0,0,0,uhat);
//      nonOp.update<0>();
//      x += nonOp.deriv().solve(-nonOp.value()),
//      nonLinearSolver.solve(x);
//      x += nonOp.value().solve(PredictRHS);
      for (int ls = 0; ls < loadSteps_; ++ls) {
        this->notify(ControlMessages::STEP_STARTED);
        nonLinearSolver.solve(x);
        this->notify(ControlMessages::SOLUTION_CHANGED);
        loadParameter += stepSize_;
        this->notify(ControlMessages::STEP_ENDED);
      }
    }

    void subscribeToNonLinearSolver(std::shared_ptr<IObserver<NonLinearSolverMessages>> observer) {
      nonLinearSolver.subscribeAll(observer);
    }

  private:
//    FEManager* feManager_;
//    DirichletManager* dirichletManager_;
    LinearAlgebraFunctions<LinearAlgebraFunctionArgs...> linearAlgebraFunctions_;
    //    FEParameterValuePair time_;
    //    using NonLinearSolver
    NonLinearSolver nonLinearSolver;
    int loadSteps_;
    double parameterBegin_;
    double parameterEnd_;
    double stepSize_;
  };

  //  template <template <typename, typename, typename> class NonLinearSolver, typename FEManager, typename ScalarType,
  //            typename DirichletManager, typename... LinearAlgebraFunctionArgs>
  //  auto makeLoadControl(FEManager& feManager, DirichletManager& dirichletManager,
  //                       const LinearAlgebraFunctions<LinearAlgebraFunctionArgs...>& linearAlgebraFunctions,
  //                       Ikarus::ILinearSolver<ScalarType>&& linearSolver, const int loadSteps,
  //                       const std::pair<double, double>& tbeginEnd) {
  //    return LoadControl<NonLinearSolver, FEManager, ScalarType, DirichletManager, LinearAlgebraFunctionArgs...>(
  //        feManager, dirichletManager, linearAlgebraFunctions, std::move(linearSolver), loadSteps, tbeginEnd);
  //  }
}  // namespace Ikarus
