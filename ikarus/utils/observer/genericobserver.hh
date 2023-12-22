// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <concepts>
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {
  template <typename Messages>
  class GenericObserver : public IObserver<Messages> {
  public:
    template <typename F>
    GenericObserver(Messages p_message, F&& p_f) : message{p_message}, f{p_f} {}
    void updateImpl(Messages p_message) override {
      if (p_message == message) {
        f(step);
        ++step;
      }
    }

    void updateImpl(Messages, double) override {}

    void updateImpl(Messages, const Eigen::VectorXd&) override {}
    Messages message;
    std::function<void(int)> f;
    int step{0};
  };
}  // namespace Ikarus
