

#pragma once

#include <dune/alugrid/grid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/nurbsgrid.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>
#include <vector>
namespace Grids {
  struct Yasp {};
  struct Alu {};
  struct Iga {};
}  // namespace Grids

template <typename GridType>
auto createGrid([[maybe_unused]] int elex = 10, [[maybe_unused]] int eley = 10) {
  //  //  /// ALUGrid Example
  if constexpr (std::is_same_v<GridType, Grids::Alu>) {
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    auto grid  = Dune::GmshReader<Grid>::read("testFiles/unstructuredTrianglesfine.msh", false);
    grid->globalRefine(0);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::Yasp>) {
    using Grid     = Dune::YaspGrid<2>;
    const double L = 1;
    const double h = 1;

    Dune::FieldVector<double, 2> bbox       = {L, h};
    std::array<int, 2> elementsPerDirection = {elex, eley};
    auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::Iga>) {
    constexpr auto dimworld        = 2;
    const std::array<int, 2> order = {2, 2};

    const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

    using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

    const std::vector<std::vector<ControlPoint>> controlPoints
        = {{{.p = {0, 0}, .w = 5}, {.p = {0.5, 0}, .w = 1}, {.p = {1, 0}, .w = 1}},
           {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 10}, {.p = {1, 0.5}, .w = 1}},
           {{.p = {0, 1}, .w = 1}, {.p = {0.5, 1}, .w = 1}, {.p = {1, 1}, .w = 1}}};

    std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

    auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
    using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

    Dune::IGA::NURBSPatchData<2, dimworld> patchData;
    patchData.knotSpans     = knotSpans;
    patchData.degree        = order;
    patchData.controlPoints = controlNet;
    auto grid               = std::make_shared<Grid>(patchData);
    grid->globalRefine(1);
    return grid;
  }
}


template <int size>
class CornerFactory {
public:
  static void construct(std::vector<Dune::FieldVector<double,size>>& values, const int corners = 10) {
    values.resize(corners);
    std::generate(values.begin(),values.end(), []() { return Ikarus::createRandomVector<Dune::FieldVector<double,size>>(); });
  }
};