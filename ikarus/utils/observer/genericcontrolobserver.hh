// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {
  class GenericControlObserver : public IObserver<ControlMessages> {
  public:
    template <typename F>
    GenericControlObserver(ControlMessages p_message, F&& p_f) : message{p_message}, f{p_f} {}
    void updateImpl(ControlMessages p_message) override {
      if (p_message == message) {
        f(step);
        ++step;
      }
    }

    void updateImpl(ControlMessages, double) override {}

    void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}
    ControlMessages message;
    std::function<void(int)> f;
    int step{0};
  };
}  // namespace Ikarus
