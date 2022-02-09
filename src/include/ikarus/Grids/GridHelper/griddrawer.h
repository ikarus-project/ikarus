//
// Created by alex on 24.06.21.
//

#pragma once

#include <matplot/matplot.h>
#include <ranges>
#include <set>
#include <dune/geometry/dimension.hh>
#include <ikarus/Grids/dunegridHelper.h>

#include "ikarus/Grids/EntityHelperFunctions.h"
#include "ikarus/Variables/VariableDefinitions.h"





template <typename GridView>
void draw(const GridView& gridView) {
  using namespace matplot;
  auto f  = figure(true);
  auto ax = gca();
  hold(ax, true);
  if constexpr (GridView::dimensionworld == 3) {
    for (auto&& edge : subentities(gridView,Dune::Dim<1>())) {
      std::array<double, 2> xEdge{}, yEdge{}, zEdge{};
      for (int i = 0; i< 2; ++i) {
        auto vertCoords = edge.geometry().corner(i);
        xEdge[i] = vertCoords[0];
        yEdge[i] = vertCoords[1];
        zEdge[i] = vertCoords[2];
      }

      auto l = ax->plot3(xEdge, yEdge, zEdge, "-o");
      l->line_width(2);
      l->color("black");
      l->marker_size(10);
      l->marker_face_color("red");
    }
  } else if constexpr (GridView::dimensionworld == 2) {
    for (auto&& edge : edges(gridView)) {
      std::array<double, 2> xEdge{}, yEdge{};
      for (int i = 0; i< 2; ++i) {
        auto vertCoords = edge.geometry().corner(i);
        xEdge[i] = vertCoords[0];
        yEdge[i] = vertCoords[1];
      }

      auto l = ax->plot(xEdge, yEdge, "-o");
      l->line_width(2);
      l->color("black");
      l->marker_size(10);
      l->marker_face_color("red");
    }
  }

  f->draw();
}

template <typename GridView, typename FEManager>
void drawDeformed(const GridView& gridView, const FEManager& feManager) {
  using namespace matplot;
  auto f1 = figure(true);
  auto ax = gca();
  hold(ax, true);
  auto fe = feManager.getFeContainer().begin();

  std::array<std::array<int, 2>, 4> edgeID{};
  edgeID[0] = {0, 2};
  edgeID[1] = {0, 1};
  edgeID[2] = {1, 3};
  edgeID[3] = {2, 3};
  //  auto eleVars    = eleVarsVec.begin();
  //  assert(surfaces(gridView).size() == feManager.size());

  auto eleVec = surfaces(gridView);
  for (auto&& ele = eleVec.begin(); ele != eleVec.end(); ++ele) {
    const auto& feVars                   = feManager.elementVariables(*fe);
    const auto& vertexDisplacementsOfEle = feVars.get(Ikarus::EntityType::vertex);
    for (int edgei = 0; auto&& edge : edges(*ele)) {
      std::array<double, 2> xEdge{}, yEdge{};
      for (int i = 0; i< 2; ++i) {
        auto vertCoords = edge.geometry().corner(i);
        xEdge[i] = vertCoords[0] + getValue(vertexDisplacementsOfEle[edgeID[edgei][i]])[0];
        yEdge[i] = vertCoords[1] + getValue(vertexDisplacementsOfEle[edgeID[edgei][i]])[1];
      }
      auto l = ax->plot(xEdge, yEdge, "-o");
      l->line_width(2);
      l->color("black");
      l->marker_size(10);
      l->marker_face_color("red");
      ++edgei;
    }
    ++fe;
  }
  f1->draw();
  //  sleep(5);
}
