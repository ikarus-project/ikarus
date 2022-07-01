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

#include <numbers>

#include <dune/alugrid/grid.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/utils/drawing/griddrawer.hh>

int main() {
  constexpr int gridDim = 2;  // (1)
  using Grid            = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
  auto grid             = Dune::GmshReader<Grid>::read("../../examples/src/testFiles/circleCoarse.msh", false);
  auto gridView         = grid->leafGridView();  // (2)

  draw(gridView);

  /// Calculate area from volume function of elements
  double area1 = 0.0;
  for (auto& element : elements(gridView))
    area1 += element.geometry().volume();

  /// Integrate function using integration rule on grid
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

  /// Naive refinement of grid and compare calculated area to pi
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
  /// write element areas to vtk
  std::vector<double> areas;
  areas.resize(gridView.size(0));

  auto& indexSet = gridView.indexSet();
  for (auto& ele : elements(gridView))
    areas[indexSet.index(ele)] = ele.geometry().volume();

  Dune::VTKWriter vtkWriter(gridView);
  vtkWriter.addCellData(areas, "area", 1);
  vtkWriter.write("TestGridEntitites");

  /// Calculate circumference and compare to pi
  double circumference = 0.0;
  for (auto& element : elements(gridView))
    if (element.hasBoundaryIntersections())
      for (auto& intersection : intersections(gridView, element))
        if (intersection.boundary()) circumference += intersection.geometry().volume();

  std::cout << circumference << " " << std::numbers::pi << std::endl;
}