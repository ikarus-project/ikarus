//
// Created by alex on 3/15/22.
//

#pragma once
#include <dune/grid/common/boundarysegment.hh>

struct UnitCircleBoundary : Dune::BoundarySegment<2, 2, double> {
  UnitCircleBoundary(const Dune::FieldVector<double, 2>& a, const Dune::FieldVector<double, 2>& b) : corners{{a, b}} {}
  Dune::FieldVector<double, 2> operator()(const Dune::FieldVector<double, 1>& local) const override {
    Dune::FieldVector<double, 2> result = {0, 0};
    double omega                        = std::acos(corners[0] * corners[1]);
    return std::sin((1 - local[0]) * omega) / sin(omega) * corners[0] + sin(local[0] * omega) / sin(omega) * corners[1];
  }

  std::array<Dune::FieldVector<double, 2>, 2> corners;
};


