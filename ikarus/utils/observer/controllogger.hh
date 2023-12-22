// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "observer.hh"
#include "observermessages.hh"

#include <chrono>

namespace Ikarus {
  class ControlLogger : public IObserver<ControlMessages> {
  public:
    void updateImpl(Ikarus::ControlMessages message) final;
    void updateImpl(Ikarus::ControlMessages, double) final {}
    void updateImpl(Ikarus::ControlMessages, int) final {}
    void updateImpl(Ikarus::ControlMessages message, const std::string& val) final;
    void updateImpl(Ikarus::ControlMessages message, int val1, const std::string& val2) final;
    void updateImpl(Ikarus::ControlMessages message, int val1, double val2) final;
    void updateImpl(Ikarus::ControlMessages, const Eigen::VectorXd&) final {}

  private:
    using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
    TimePoint start;
    TimePoint stop;
    std::chrono::milliseconds duration;
  };
}  // namespace Ikarus