//
// Created by lex on 17/12/2021.
//

#pragma once
#include <memory>

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/Variables/ParameterFactory.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include <ikarus/utils/Observer/observer.h>
#include <ikarus/utils/Observer/observerMessages.h>

namespace Ikarus {

  template <typename FEManager, typename DirichletManager, typename... LinearAlgebraFunctionArgs>
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(FEManager& feManager,DirichletManager& dirichletManager,
                const LinearAlgebraFunctions<LinearAlgebraFunctionArgs...>& linearAlgebraFunctions, const int loadSteps,
                const std::pair<double, double>& tbeginEnd)
        : feManager_{&feManager},
          dirichletManager_{&dirichletManager},
          linearAlgebraFunctions_{linearAlgebraFunctions},
          loadSteps_{loadSteps},
          tBegin_{tbeginEnd.first},
          tEnd_{tbeginEnd.second},
          timeStepSize_{(tEnd_ - tBegin_) / loadSteps_} {}

    void newtonRaphson( auto& nonLinearOperator, auto& x)
    {
      const double tol   = 1e-10;
      const int maxInter = 2;
      auto rNorm = nonLinearOperator.value().norm();
      int iter   = 0;
      while (rNorm > tol && iter <= maxInter) {
        const auto& residual = nonLinearOperator.value();
        const auto& K        = nonLinearOperator.derivative();

        const Eigen::VectorXd D = -K.ldlt().solve(residual);
        x += dirichletManager_->viewAsFullVector(D);
        nonLinearOperator.updateAll();
        rNorm = residual.norm();
        this->notify(ControlMessages::SOLUTION_CHANGED);
        ++iter;
      }
    }

    void run() {
      this->notify(ControlMessages::CONTROL_STARTED);
      auto time = Ikarus::FEParameterFactory::createParameter(Ikarus::FEParameter::time, 1);
      time.value[0] = timeStepSize_;
      Ikarus::NonLinearOperator nonLinearOperator(linearAlgebraFunctions_, parameter(time));
      auto& x            = feManager_->getVariables();

      NewtonRaphson nr(nonLinearOperator,[this](auto&&x, auto&& D){ x+=dirichletManager_->viewAsFullVector(D);});
      for (int ls = 0; ls < loadSteps_; ++ls) {
        nr.solve(x);
        std::cout<<x<<std::endl;
        time.value[0] += timeStepSize_;
        this->notify(ControlMessages::LOADSTEP_ENDED);
      }
    }

  private:
    FEManager* feManager_;
    DirichletManager* dirichletManager_;
    LinearAlgebraFunctions<LinearAlgebraFunctionArgs...> linearAlgebraFunctions_;
    int loadSteps_;
    double tBegin_;
    double tEnd_;
    double timeStepSize_;
  };
}  // namespace Ikarus
