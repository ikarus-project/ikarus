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

  //  template <typename... Args>
  //  struct NonLinearSolverArguments {
  //    std::tuple<std::reference_wrapper<std::remove_cvref_t<Args>>...> args;
  //  };
  //
  //  template <typename... Args>
  //  auto nonlinearSolverArguments(Args&&... args) {
  //    return NonLinearSolverArguments<Args&&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
  //  }

  //  template <template <typename, typename, typename> class NonLinearSolver, typename FEManager,
  //            typename DirichletManager, typename TypeListOne, typename TypeListTwo>
  //  class LoadControl {
  //  public:
  //    LoadControl(FEManager& feManager, DirichletManager& dirichletManager, const TypeListOne&
  //    nonLinearSolverArguments,
  //                const TypeListTwo& linearAlgebraFunctions, const int loadSteps,
  //                const std::pair<double, double>& tbeginEnd) {}
  //  };

  template <typename NonLinearSolver, typename FEManager, typename DirichletManager,
            typename... LinearAlgebraFunctionArgs>
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(FEManager& feManager, DirichletManager& dirichletManager, NonLinearSolver&& p_nonLinearSolver,
                const int loadSteps, const std::pair<double, double>& tbeginEnd)
        : feManager_{&feManager},
          dirichletManager_{&dirichletManager},
          nonLinearSolver{std::move(p_nonLinearSolver)},
          loadSteps_{loadSteps},
          parameterBegin_{tbeginEnd.first},
          parameterEnd_{tbeginEnd.second},
          stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {}

    void run() {
      this->notify(ControlMessages::CONTROL_STARTED);
      auto& loadParameter    = nonLinearSolver.nonLinearOperator().template nthParameter<0>();
      loadParameter.value[0] = 0.0;
      auto& x                = feManager_->getVariables();
      for (int ls = 0; ls < loadSteps_ + 1; ++ls) {
        nonLinearSolver.solve(x);
        this->notify(ControlMessages::SOLUTION_CHANGED);
        loadParameter.value[0] += stepSize_;
        this->notify(ControlMessages::LOADSTEP_ENDED);
      }
    }

    void subscribeToNonLinearSolver(std::shared_ptr<IObserver<NonLinearSolverMessages>> observer) {
      nonLinearSolver.subscribeAll(observer);
    }

  private:
    FEManager* feManager_;
    DirichletManager* dirichletManager_;
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
