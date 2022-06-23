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

#include <dune/geometry/dimension.hh>

template <typename GridView>
void draw(const GridView& gridView) {
  using namespace matplot;
  auto f  = figure(true);
  auto ax = gca();
  hold(ax, true);
  if constexpr (GridView::dimensionworld == 3) {
    for (auto&& edge : subentities(gridView, Dune::Dim<1>())) {
      std::array<double, 2> xEdge{}, yEdge{}, zEdge{};
      for (int i = 0; i < 2; ++i) {
        auto vertCoords = edge.geometry().corner(i);
        xEdge[i]        = vertCoords[0];
        yEdge[i]        = vertCoords[1];
        zEdge[i]        = vertCoords[2];
      }

      auto l = ax->plot3(xEdge, yEdge, zEdge, "-o");
      l->line_width(2);
      l->color("black");
      l->marker_size(10);
      l->marker_face_color("red");
    }
  } else if constexpr (GridView::dimensionworld == 2) {  // FIXME reduce code duplciation
    for (auto&& edge : edges(gridView)) {
      std::array<double, 2> xEdge{}, yEdge{};
      for (int i = 0; i < 2; ++i) {
        auto vertCoords = edge.geometry().corner(i);
        xEdge[i]        = vertCoords[0];
        yEdge[i]        = vertCoords[1];
      }

      auto l = ax->plot(xEdge, yEdge, "-o");
      l->line_width(2);
      l->color("black");
      l->marker_size(10);
      l->marker_face_color("red");
    }
  }

  else if constexpr (GridView::dimensionworld == 1) {
    for (auto&& edge : elements(gridView)) {
      std::array<double, 2> xEdge{}, yEdge{};
      for (int i = 0; i < 2; ++i) {
        auto vertCoords = edge.geometry().corner(i);
        xEdge[i]        = vertCoords[0];
        yEdge[i]        = 0.0;
      }

      auto l = ax->plot(xEdge, yEdge, "-o");
      l->line_width(2);
      l->color("black");
      l->marker_size(10);
      l->marker_face_color("red");
    }
  }

  f->show();
  f.reset();
}
