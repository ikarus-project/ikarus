// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/Observer/observer.hh>
#include <ikarus/utils/Observer/observermessages.hh>
#include <ikarus/utils/drawing/griddrawer.hh>

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

  void updateImpl(ControlMessages, double) override {}

  void updateImpl(ControlMessages, const Eigen::VectorXd&) override {}

private:
  GridView const* gridView_;
  FEManager const* feManager_;
};
