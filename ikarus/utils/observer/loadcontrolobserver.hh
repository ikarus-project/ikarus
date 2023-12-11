// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <spdlog/spdlog.h>

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {
  class LoadControlObserver : public IObserver<ControlMessages> {
  public:
    void updateImpl(ControlMessages message) override {
      switch (message) {
        case ControlMessages::CONTROL_STARTED: {
          spdlog::info("Load control ======= STARTED  =======");
        } break;
        case ControlMessages::CONTROL_ENDED: {
          spdlog::info("Load control ======= ENDED  =======");
        } break;
        default:
          break;  //   default: do nothing when notified
      }
    }

    using IObserver::updateImpl;
    void updateImpl(ControlMessages, double) override {}
    void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}
  };
}  // namespace Ikarus
