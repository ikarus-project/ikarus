//
// Created by lex on 14/12/2021.
//

#pragma once
#include <string>

#include "spdlog/spdlog.h"

#include "ikarus/Grids/GridHelper/griddrawer.h"
#include "ikarus/utils/Observer/observer.h"
#include "ikarus/utils/Observer/observerMessages.h"

template <typename GridView, typename FEManager>
class GridDrawerObserver : public IObserver<ControlMessages> {
public:
  GridDrawerObserver(const GridView& gridView, const FEManager& feManager)
      : gridView_{&gridView}, feManager_{&feManager} {}

  void updateImpl(ControlMessages message) override {
    switch (message) {
      case ControlMessages::SOLUTION_CHANGED:
        drawDeformed(*gridView_, *feManager_);
        break;
      default:
        break;  //   default: do nothing when notified
    }
  }

  void updateImpl(ControlMessages message, double val) override {}

  void updateImpl(ControlMessages message, const Eigen::VectorXd& vec) override {}

private:
  GridView const* gridView_;
  FEManager const* feManager_;
};
