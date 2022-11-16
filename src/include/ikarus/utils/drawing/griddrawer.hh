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

#include <matplot/matplot.h>
#include <ranges>
#include <set>
#include <thread>

#include <dune/geometry/dimension.hh>

template <typename GridView>
void draw(const GridView& gridView, bool forever = false) {
  using namespace matplot;
  auto f  = figure(true);
  auto ax = gca();
  hold(ax, true);
  constexpr int edgeCodim = GridView::dimension - 1;
  for (auto&& element : elements(gridView)) {
    std::array<std::array<double, 2>, GridView::dimensionworld> edgeCoords{};
    for (size_t edgeIndex = 0; edgeIndex < element.subEntities(edgeCodim); ++edgeIndex) {
      auto edge = element.template subEntity<edgeCodim>(edgeIndex);
      for (int i = 0; i < 2; ++i) {
        const auto vertCoords = edge.geometry().corner(i);
        for (int j = 0; j < GridView::dimensionworld; ++j)
          edgeCoords[j][i] = vertCoords[j];
      }
      if constexpr (GridView::dimensionworld == 3) {
        auto l = ax->plot3(edgeCoords[0], edgeCoords[1], edgeCoords[2], "-o");
        l->line_width(2);
        l->color("black");
        l->marker_size(10);
        l->marker_face_color("red");
      } else {
        auto l = ax->plot(edgeCoords[0], edgeCoords[1], "-o");
        l->line_width(2);
        l->color("black");
        l->marker_size(10);
        l->marker_face_color("red");
      }
    }
  }

  if (forever)
    f->show();
  else {
    f->draw();
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(5s);
  }
  f.reset();
}