//
// Created by Alex on 21.07.2021.
//

#include <numbers>

#include <dune/alugrid/grid.hh>
#include <dune/geometry/quadraturerules.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/Grids/GridHelper/griddrawer.h>

int main() {
  constexpr int gridDim = 2;
  using Grid            = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
  auto grid             = Dune::GmshReader<Grid>::read("../../examples/src/testFiles/circleCoarse.msh", false);
  auto gridView         = grid->leafGridView();

  draw(gridView);

  double area1 = 0.0;
  for (auto& element : elements(gridView)) {
    area1 += element.geometry().volume();
  }

  auto f       = [](auto&& global) { return sqrt(global[0] * global[0] + global[1] * global[1]); };
  double area2 = 0.0;
  for (auto& element : elements(gridView)) {
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(element.type(), 1, Dune::QuadratureType::GaussLegendre);
    for (auto& gp : rule)
      //      area2 += element.geometry().integrationElement(gp.position()) * gp.weight();
      area2 += f(element.geometry().global(gp.position())) * element.geometry().integrationElement(gp.position())
               * gp.weight();  // integrationElement --> JacobiDeterminant
  }

  std::cout << area1 << " " << area2 << std::endl;

  for (int i = 0; i < 3; ++i) {
    area1 = 0.0;
    grid->globalRefine(1);

    auto gridViewRefined = grid->leafGridView();
    std::cout << "This gridview contains: ";
    std::cout << gridViewRefined.size(0) << " elements" << std::endl;
    draw(gridViewRefined);
    for (auto& element : elements(gridViewRefined)) {
      area1 += element.geometry().volume();
    }
    std::cout << area1 << " " << std::numbers::pi << std::endl;
  }

  double circumference = 0.0;
  for (auto& element : elements(gridView))
    if (element.hasBoundaryIntersections())
      for (auto& intersection : intersections(gridView, element))
        if (intersection.boundary()) circumference += intersection.geometry().volume();

  std::cout << circumference << std::endl;
}