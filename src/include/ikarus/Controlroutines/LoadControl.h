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

  template <template <typename, typename, typename> class NonLinearSolver, typename FEManager, typename ScalarType,
            typename DirichletManager, typename... LinearAlgebraFunctionArgs>
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(FEManager& feManager, DirichletManager& dirichletManager,
                const LinearAlgebraFunctions<LinearAlgebraFunctionArgs...>& linearAlgebraFunctions,
                Ikarus::ILinearSolver<ScalarType>&& linearSolver, const int loadSteps,
                const std::pair<double, double>& tbeginEnd)
        : feManager_{&feManager},
          dirichletManager_{&dirichletManager},
          linearAlgebraFunctions_{linearAlgebraFunctions},
          time_{Ikarus::FEParameterFactory::createParameter(Ikarus::FEParameter::time, 1)},
          nonLinearSolver{NonLinearOperatorType(linearAlgebraFunctions_, parameter(time_)), std::move(linearSolver),
                          [this](decltype(feManager_->getVariables())& x, const Eigen::VectorX<ScalarType>& D) {
                            x += dirichletManager_->viewAsFullVector(D);
                          }},
          loadSteps_{loadSteps},
          tBegin_{tbeginEnd.first},
          tEnd_{tbeginEnd.second},
          timeStepSize_{(tEnd_ - tBegin_) / loadSteps_} {}

    void run() {
      this->notify(ControlMessages::CONTROL_STARTED);
      time_.value[0] = timeStepSize_;
      auto& x        = feManager_->getVariables();
      for (int ls = 0; ls < loadSteps_; ++ls) {
        nonLinearSolver.solve(x);
        time_.value[0] += timeStepSize_;
        this->notify(ControlMessages::LOADSTEP_ENDED);
      }
    }

  private:
    FEManager* feManager_;
    DirichletManager* dirichletManager_;
    LinearAlgebraFunctions<LinearAlgebraFunctionArgs...> linearAlgebraFunctions_;
    FEParameterValuePair time_;
    using NonLinearOperatorType
        = Ikarus::NonLinearOperator<LinearAlgebraFunctions<LinearAlgebraFunctionArgs...>, decltype(parameter(time_))>;
    NonLinearSolver<NonLinearOperatorType, Ikarus::ILinearSolver<ScalarType>,
                    std::function<void(decltype(feManager_->getVariables())&, const Eigen::VectorX<ScalarType>&)>>
        nonLinearSolver;
    int loadSteps_;
    double tBegin_;
    double tEnd_;
    double timeStepSize_;
  };

  template <template <typename, typename, typename> class NonLinearSolver, typename FEManager, typename ScalarType,
            typename DirichletManager, typename... LinearAlgebraFunctionArgs>
  auto makeLoadControl(FEManager& feManager, DirichletManager& dirichletManager,
                       const LinearAlgebraFunctions<LinearAlgebraFunctionArgs...>& linearAlgebraFunctions,
                       Ikarus::ILinearSolver<ScalarType>&& linearSolver, const int loadSteps,
                       const std::pair<double, double>& tbeginEnd) {
    return LoadControl<NonLinearSolver, FEManager, ScalarType, DirichletManager, LinearAlgebraFunctionArgs...>(
        feManager, dirichletManager, linearAlgebraFunctions, std::move(linearSolver), loadSteps, tbeginEnd);
  }
}  // namespace Ikarus
