//
// Created by lex on 14/12/2021.
//

#pragma once
#include <string>

#include "spdlog/spdlog.h"

#include "ikarus/utils/Observer/observer.h"
#include "ikarus/utils/Observer/observerMessages.h"



class ControlLogger : public IObserver<ControlMessages> {
public:
  void updateImpl(ControlMessages message) override {
    switch (message) {
      case ControlMessages::CONTROL_STARTED:
        spdlog::info("Control started");
        break;
      case ControlMessages::ITERATION_ENDED:
        spdlog::info("Iteration has ended");
        break;
      case ControlMessages::LOADSTEP_ENDED:
        spdlog::info("============================================\n");
        spdlog::info("Loadstep has ended");
        spdlog::info("============================================\n");
        break;
      case ControlMessages::SOLUTION_CHANGED:
        spdlog::info("ControlMessages::SOLUTION_CHANGED");
        break;
      default:
        break;  //   default: do nothing when notified
    }
  }

  void updateImpl(ControlMessages message, double val) override {
    switch (message) {
      case ControlMessages::RESIDUALNORM_UPDATED:
        spdlog::info("Residual norm is {:03.2f}", val);
        break;
      default:
        break;
    }
  }

  void updateImpl(ControlMessages message, const Eigen::VectorXd& vec) override {

  }
};
