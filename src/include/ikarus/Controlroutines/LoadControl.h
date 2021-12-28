//
// Created by lex on 17/12/2021.
//

#pragma once
#include <memory>

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/Variables/ParameterFactory.h"
#include <ikarus/utils/Observer/observer.h>
#include <ikarus/utils/Observer/observerMessages.h>

namespace Ikarus {

  template <typename FEManager, typename... LinearAlgebraFunctionArgs>
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(FEManager& feManager,
                const LinearAlgebraFunctions<LinearAlgebraFunctionArgs...>& linearAlgebraFunctions, const int loadSteps,
                const std::pair<double, double>& tbeginEnd)
        : feManager_{&feManager},
          linearAlgebraFunctions_{linearAlgebraFunctions},
          loadSteps_{loadSteps},
          tBegin_{tbeginEnd.first},
          tEnd_{tbeginEnd.second},
          timeStepSize_{(tEnd_ - tBegin_) / loadSteps_} {}

//    template<typename LoadFunction>
//    void setLoadFunction(LoadFunction&& loadfunction)
//                         {
//      loadFunction_         = loadfunction;
//
//}

    void run() {
      this->notify(ControlMessages::CONTROL_STARTED);
      auto lambda = Ikarus::FEParameterFactory::createParameter(Ikarus::FEParameter::time, 1);
      Eigen::VectorXd lambdav(1);
      lambdav[0] = 0.1;
      setValue(lambda.value, lambdav);
      Ikarus::NonLinearOperator nonLinearOperator(linearAlgebraFunctions_, parameter(lambda));
      auto& x            = feManager_->getVariables();
      const double tol   = 1e-10;
      const int maxInter = 2;
      for (int ls = 0; ls < loadSteps_; ++ls) {
        auto rNorm = nonLinearOperator.value().norm();
        int iter   = 0;
        while (rNorm > tol && iter <= maxInter) {
          const auto& residual = nonLinearOperator.value();
          const auto& K        = nonLinearOperator.derivative();

          const Eigen::VectorXd D = -K.ldlt().solve(residual);
          x += D;
          nonLinearOperator.updateAll();
          rNorm = residual.norm();
          this->notify(ControlMessages::RESIDUALNORM_UPDATED, rNorm);
          this->notify(ControlMessages::SOLUTION_CHANGED);
          this->notify(ControlMessages::ITERATION_ENDED);
          ++iter;
        }
        lambdav[0] += 0.1;
        setValue(lambda.value, lambdav);
        this->notify(ControlMessages::LOADSTEP_ENDED);
        this->notify(ControlMessages::LOADSTEP_ENDED);
      }
    }

  private:
    FEManager* feManager_;
    LinearAlgebraFunctions<LinearAlgebraFunctionArgs...> linearAlgebraFunctions_;
    int loadSteps_;
    double tBegin_;
    double tEnd_;
    double timeStepSize_;
  };
}  // namespace Ikarus
