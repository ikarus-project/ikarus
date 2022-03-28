//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/alugrid/grid.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/boundarysegment.hh>
//#include <dune/common/function.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/utils/drawing/griddrawer.h"

struct UnitCircleBoundary : Dune::BoundarySegment<2, 2, double> {
  UnitCircleBoundary(const Dune::FieldVector<double, 2>& a, const Dune::FieldVector<double, 2>& b) : corners{{a, b}} {}
  Dune::FieldVector<double, 2> operator()(const Dune::FieldVector<double, 1>& local) const override {
    Dune::FieldVector<double, 2> result = {0, 0};
    double omega                        = std::acos(corners[0] * corners[1]);
    return std::sin((1 - local[0]) * omega) / sin(omega) * corners[0] + sin(local[0] * omega) / sin(omega) * corners[1];
  }

  std::array<Dune::FieldVector<double, 2>, 2> corners;
};

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);

  /// Create grid from 6 triagnles align in unit disc
  using namespace Dune;
  constexpr int gridDim = 2;
  Dune::GridFactory<Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>> gridFactory;
  //  std::array<FieldVector<double, 2>, 4> corners0 = {{{-sqrt(2) / 2, -sqrt(2) / 2}, {sqrt(2) / 2, -sqrt(2) / 2},
  //  {sqrt(2) / 2, sqrt(2) / 2}, {-sqrt(2) / 2, sqrt(2) / 2}}};
  Eigen::Vector2d v(1, 0);
  std::array<FieldVector<double, 2>, 6> corners0;
  Eigen::Rotation2D<double> R;
  R.angle() = 0.0;
  for (auto& corner : corners0) {
    Eigen::Vector2d a = R * v;
    corner[0]         = a[0];
    corner[1]         = a[1];
    R.angle() += 60.0 / 180.0 * std::numbers::pi;
  }

  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex(corners0[0]);
  gridFactory.insertVertex(corners0[1]);
  gridFactory.insertVertex(corners0[2]);
  gridFactory.insertVertex(corners0[3]);
  gridFactory.insertVertex(corners0[4]);
  gridFactory.insertVertex(corners0[5]);

  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 1, 2});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 2, 3});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 3, 4});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 4, 5});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 5, 6});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 6, 1});

  /// Create boundary segments which map the boundaries onto the unit circle
  gridFactory.insertBoundarySegment({1, 2}, std::make_shared<UnitCircleBoundary>(corners0[0], corners0[1]));
  gridFactory.insertBoundarySegment({2, 3}, std::make_shared<UnitCircleBoundary>(corners0[1], corners0[2]));
  gridFactory.insertBoundarySegment({3, 4}, std::make_shared<UnitCircleBoundary>(corners0[2], corners0[3]));
  gridFactory.insertBoundarySegment({4, 5}, std::make_shared<UnitCircleBoundary>(corners0[3], corners0[4]));
  gridFactory.insertBoundarySegment({5, 6}, std::make_shared<UnitCircleBoundary>(corners0[4], corners0[5]));
  gridFactory.insertBoundarySegment({6, 1}, std::make_shared<UnitCircleBoundary>(corners0[5], corners0[0]));

  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();
  draw(gridView);

  double area1 = 0.0;
  for (const auto& element : elements(gridView)) {
    area1 += element.geometry().volume();
  }

  auto f       = [](auto&& global) { return sqrt(global[0] * global[0] + global[1] * global[1]); };
  double area2 = 0.0;
  for (auto& element : elements(gridView)) {
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(element.type(), 1, Dune::QuadratureType::GaussLegendre);
    for (auto& gp : rule)
      area2 += element.geometry().integrationElement(gp.position()) * gp.weight();
  }

  std::cout << "Area2 " << area2 << " " << std::numbers::pi << std::endl;
  for (int i = 0; i < 10; ++i) {
    area1 = 0.0;
    /// Refine grid entities if they live at the boundary
    //    grid->globalRefine(1);
    for (const auto& ele : elements(grid->leafGridView())) {
      if (ele.hasBoundaryIntersections()) grid->mark(1, ele);
    }
    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();
    auto gridViewRefined = grid->leafGridView();

    std::cout << "This gridview contains: ";
    std::cout << gridViewRefined.size(0) << " elements" << std::endl;

    for (auto& element : elements(gridViewRefined))
      area1 += element.geometry().volume();

    std::cout << "area1 " << area1 << " " << std::numbers::pi << std::endl;
    draw(gridViewRefined);
  }
  /// Calculate circumference and compare to pi
  double circumference = 0.0;
  for (auto& element : elements(grid->leafGridView()))
    if (element.hasBoundaryIntersections())
      for (auto& intersection : intersections(grid->leafGridView(), element))
        if (intersection.boundary()) circumference += intersection.geometry().volume();

  std::cout << circumference / 2 << " " << std::numbers::pi << std::endl;
}