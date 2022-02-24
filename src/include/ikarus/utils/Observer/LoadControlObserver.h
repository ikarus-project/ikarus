//
// Created by lex on 14/12/2021.
//

#pragma once
#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "spdlog/spdlog.h"

#include "ikarus/utils/Observer/observer.h"
#include "ikarus/utils/Observer/observerMessages.h"

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

    void updateImpl(ControlMessages, double) override {}
    void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}
  };
}  // namespace Ikarus