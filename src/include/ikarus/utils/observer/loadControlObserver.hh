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
#include <string>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <spdlog/spdlog.h>

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observerMessages.hh>

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