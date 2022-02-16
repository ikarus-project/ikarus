//
// Created by lex on 14/12/2021.
//

#pragma once
#include <concepts>
#include <string>

#include "spdlog/spdlog.h"

#include "ikarus/utils/Observer/observer.h"
#include "ikarus/utils/Observer/observerMessages.h"

namespace Ikarus {
  class GenericControlObserver : public IObserver<ControlMessages> {
  public:
    template<typename F>
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