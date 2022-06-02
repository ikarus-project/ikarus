//
// Created by alex on 24.06.21.
//

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
    for (auto&& ele : elements(gridView)) {
      std::vector<double> xEdge( ele.geometry().corners()), yEdge( ele.geometry().corners());
      for (int i = 0; i < ele.geometry().corners(); ++i) {
        auto vertCoords = ele.geometry().corner(i);
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
