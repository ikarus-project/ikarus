// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/Observer/observer.hh>
#include <ikarus/utils/Observer/observerMessages.hh>

class ControlLogger : public IObserver<ControlMessages> {
public:
  void updateImpl(ControlMessages message) override {
    switch (message) {
      case ControlMessages::CONTROL_STARTED:
        spdlog::info("Control started");
        break;
      case ControlMessages::STEP_STARTED:
        spdlog::info("============================================");
        spdlog::info("Control step has started");
        spdlog::info("============================================");
        break;
      case ControlMessages::STEP_ENDED:
        spdlog::info("============================================");
        spdlog::info("Control step has ended");
        break;
      case ControlMessages::CONTROL_ENDED:
        spdlog::info("Control ended");
        spdlog::info("============================================");
        break;
      default:
        break;  //   default: do nothing when notified
    }
  }

  void updateImpl(ControlMessages message, double val) override {
    switch (message) {
      case ControlMessages::STEP_ENDED:
        spdlog::info("Step size is {:04.3e}", val);
        break;
      default:
        break;
    }
  }

  void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}
};
