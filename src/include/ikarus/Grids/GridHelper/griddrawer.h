//
// Created by alex on 24.06.21.
//

#pragma once

#include <matplot/matplot.h>
#include <ranges>
#include <set>
template <typename GridView>
void draw(GridView& gridView) {
  using namespace matplot;
  auto f = figure(true);
  auto ax = gca();
  hold(ax, true);
  for (auto&& edge : edges(gridView)) {
    std::array<double, 2> xEdge{}, yEdge{}, zEdge{};
    for (int i = 0; auto&& vert : vertices(edge)) {
      xEdge[i] = vert->getPosition()[0];
      yEdge[i] = vert->getPosition()[1];
      zEdge[i] = vert->getPosition()[2];
      ++i;
    }

    auto l = ax->plot3(xEdge, yEdge, zEdge, "-o");
    l->line_width(2);
    l->color("black");
    l->marker_size(10);
    l->marker_face_color("red");
  }
  f->show();
}
