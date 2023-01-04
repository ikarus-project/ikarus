// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <clipper2/clipper.h>


int main(int argc, char** argv) {
  using namespace Clipper2Lib;
  Dune::MPIHelper::instance(argc, argv);

  //set up the subject and clip polygons ...
  Paths64 subject;
  subject.push_back(Ellipse(Rect64(100,100,300,300)));
  subject.push_back(Ellipse(Rect64(125,130,275,180)));
  subject.push_back(Ellipse(Rect64(125,220,275,270)));

  Paths64 clip;
  clip.push_back(Ellipse(Rect64(140,70,220,320)));

  //get the intersection
  Paths64 solution = Intersect(subject, clip, FillRule::EvenOdd);

  //DrawPolygons
}
