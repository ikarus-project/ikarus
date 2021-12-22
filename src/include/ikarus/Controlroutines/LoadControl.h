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
                const LinearAlgebraFunctions<LinearAlgebraFunctionArgs...>& linearAlgebraFunctions, const int loadSteps)
        : feManager_{&feManager}, linearAlgebraFunctions_{linearAlgebraFunctions}, loadSteps_{loadSteps} {}

    void run() {
      this->notify(ControlMessages::CONTROL_STARTED);
      auto lambda = Ikarus::FEParameterFactory::createParameter(Ikarus::FEParameter::loadfactor, 1);
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
          auto& residual = nonLinearOperator.value();
          //          std::cout<<"residual"<<std::endl;
          //          std::cout<<residual<<std::endl;
          //          std::cout<<"residual"<<std::endl;
          auto& K = nonLinearOperator.derivative();
          for (int i = 1; i < K.cols(); ++i) {
            K.coeffRef(i, 0) = 0.0;
            K.coeffRef(i, 1) = 0.0;
            K.coeffRef(i, 3) = 0.0;
            K.coeffRef(0, i) = 0.0;
            K.coeffRef(1, i) = 0.0;
            K.coeffRef(3, i) = 0.0;
          }
          K.coeffRef(0, 0) = 1;
          K.coeffRef(1, 1) = 1;
          K.coeffRef(3, 3) = 1;

          const Eigen::VectorXd D = K.ldlt().solve(-residual);
          x += D;
          //          std::cout<<"D.transpose()"<<std::endl;
          //          std::cout<<D.transpose()<<std::endl;
          //          std::cout<<residual.transpose()<<std::endl;
          nonLinearOperator.updateAll();
          //          std::cout<<residual.transpose()<<std::endl;
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
  };
}  // namespace Ikarus
