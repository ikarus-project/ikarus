// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/Observer/observer.hh>
#include <ikarus/utils/Observer/observermessages.hh>

class ControlLogger : public IObserver<ControlMessages> {
public:
  void updateImpl(ControlMessages message) override {
    switch (message) {
      case ControlMessages::CONTROL_STARTED:
        spdlog::info("Control started");
        break;
      case ControlMessages::STEP_STARTED:
        spdlog::info("============================================");
        spdlog::info("Controlstep has started");
        spdlog::info("============================================");
        break;
      case ControlMessages::STEP_ENDED:
        spdlog::info("============================================");
        spdlog::info("Controlstep has ended");
        break;
      case ControlMessages::CONTROL_ENDED:
        spdlog::info("Control ended");
        spdlog::info("============================================");
        break;
      default:
        break;  //   default: do nothing when notified
    }
  }

  void updateImpl(ControlMessages, double) override {}

  void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}
};
