//
// Created by lex on 14/12/2021.
//

#pragma once
#include <ikarus/utils/observer.h>
#include <string>
#include "spdlog/spdlog.h"

enum class ControlMessages { BEGIN, ITERATION_ENDED, LOADSTEP_ENDED, RESIDUALNORM_UPDATED, END };



class ControlLogger : public IObserver<ControlMessages> {
public:

  void updateImpl(ControlMessages message) override {
    switch (message) {
      case ControlMessages::ITERATION_ENDED:
        spdlog::info("Iteration has ended");
        break;
      case ControlMessages::LOADSTEP_ENDED:
        spdlog::info("============================================\n");
        spdlog::info("Loadstep has ended" );
        spdlog::info("============================================\n");
        break;
      default:
        break; //   default: do nothing when notified
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

};
