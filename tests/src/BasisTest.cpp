//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <Eigen/Core>
#include <dune/grid/yaspgrid.hh>
#include <ikarus/FEManager/DefaultFEManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

TEST(FEManager, FEManagertest) {
    using Grid = Dune::YaspGrid<2>;
    const double L    = 1;
    const double h    = 1;
    const size_t elex = 1;
    const size_t eley = 1;

    Dune::FieldVector<double, 2> bbox = {L, h};
    std::array<int, 2> eles           = {elex, eley};
    auto grid                         = std::make_shared<Grid>(bbox, eles);
  auto gridView = grid->leafGridView();

  std::vector<Ikarus::FiniteElements::IFiniteElement> feContainer;

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));
  auto localView = basis.localView();
  for (auto& e :elements(gridView)) {
    localView.bind(e);

    std::cout<<localView.size()<<std::endl;
  }

}
