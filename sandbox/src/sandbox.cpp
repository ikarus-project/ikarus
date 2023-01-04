// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <clipper2/clipper.h>
#include <iostream>
#include <matplot/matplot.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/iga/nurbsgrid.hh>

#include <ikarus/utils/drawing/griddrawer.hh>

void drawPaths(auto paths, std::string&& lineColor, bool isLoop = true) {
  std::vector<double> x;
  std::vector<double> y;
  using namespace matplot;
  auto f  = gcf();
  auto ax = gca();
  hold(ax, true);
  for (auto& path : paths) {
    x.clear();
    y.clear();
    for (auto& point : path) {
      x.push_back(point.x);
      y.push_back(point.y);
    }
    if (isLoop)  // append first point to end again to close the loop
    {
      x.push_back(path.front().x);
      y.push_back(path.front().y);
    }

    auto l = ax->plot(x, y);
    l->line_width(2);
    l->color(lineColor);
    //        l->marker_size(10);
    //        l->marker_face_color(lineColor);
  }
}

int main(int argc, char** argv) {
  using namespace Clipper2Lib;
  Dune::MPIHelper::instance(argc, argv);

  PathsD subject;
  subject.push_back(Ellipse(RectD(100, 100, 300, 300), 100));
  subject.push_back(Ellipse(RectD(125, 130, 275, 180), 100));
  subject.push_back(Ellipse(RectD(125, 220, 275, 270), 100));

  PathsD clip;
  clip.push_back(Ellipse(RectD(140, 70, 220, 320), 100));
  PathsD solution = Intersect(subject, clip, FillRule::EvenOdd, 5);

  auto f = matplot::figure(true);
  drawPaths(subject, "red");
  drawPaths(clip, "black");
  drawPaths(solution, "blue");
  f->show();

  constexpr int griddim                                    = 2;
  constexpr int dimworld                                   = 2;
  const std::array<std::vector<double>, griddim> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointType;
  const double Lx    = 1;
  const double Ly    = 1;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {0, Ly}, .w = 1}}, {{.p = {Lx, 0}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, griddim> dimsize = {2, 2};

  std::vector<double> dofsVec;
  std::vector<double> l2Evcector;
  auto controlNet = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<griddim, dimworld>;
  Dune::IGA::NURBSPatchData<griddim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;
  /// Increate polynomial degree in each direction
  patchData = Dune::IGA::degreeElevate(patchData, 0, 1);
  patchData = Dune::IGA::degreeElevate(patchData, 1, 1);
  Grid grid(patchData);
  grid.globalRefine(4);
  auto gridView = grid.leafGridView();
  PathsD edges(1);
  PathsD alledges;
  PathsD clip2;
  PathsD trimmedandEmptyElements;
  clip2.push_back(Ellipse(RectD(0.1, 0.8, 0.7, 0.1), 100));

  for (int edgeIndex = 0; auto& element : Dune::elements(gridView)) {
    auto geo  = element.geometry();
    auto pos0 = geo.corner(0);
    auto pos1 = geo.corner(1);
    auto pos2 = geo.corner(3);  // see dune book page 127 Figure 5.12
    auto pos3 = geo.corner(2);
    edges.front().clear();
    edges.front().emplace_back(pos0[0], pos0[1]);
    edges.front().emplace_back(pos1[0], pos1[1]);
    edges.front().emplace_back(pos2[0], pos2[1]); //swap order of 2 and 3
    edges.front().emplace_back(pos3[0], pos3[1]);
    alledges.insert(alledges.end(), edges.begin(), edges.end());
    Clipper2Lib::ClipperD clipper(5);
    Clipper2Lib::PathsD clippedEdges;

    clipper.AddSubject(edges);
    clipper.AddClip(clip2);
    clipper.Execute(ClipType::Intersection, FillRule::NonZero, clippedEdges);
    trimmedandEmptyElements.insert(trimmedandEmptyElements.end(), clippedEdges.begin(), clippedEdges.end());
  }
  auto f2 = matplot::figure(true);
  drawPaths(alledges, "black");

  drawPaths(clip2, "red");
  drawPaths(trimmedandEmptyElements, "blue");
  f2->show();
}
