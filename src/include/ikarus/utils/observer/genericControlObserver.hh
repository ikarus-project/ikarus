/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */



#pragma once
#include <concepts>
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observerMessages.hh>

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