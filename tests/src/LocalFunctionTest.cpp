//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <array>
#include <complex>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>

#include <Eigen/Core>

#include "ikarus/Geometries/GeometryWithExternalInput.h"
#include "ikarus/LocalBasis/localBasis.h"
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/Geometries/SimpleLocalFunction.h>

const double tol = 1e-15;
using namespace Dune::Functions::BasisFactory;
TEST(LocalFunction, SimpleLocalFunction) {
  using Grid        = Dune::YaspGrid<2>;
  const double Lx   = 7;
  const double Ly   = 5;
  const size_t elex = 1;
  const size_t eley = 1;

  Dune::FieldVector<double, 2> bbox = {Lx, Ly};
  std::array<int, 2> eles           = {elex, eley};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  auto gridView      = grid->leafGridView();
  auto basis = makeBasis(gridView, power<3>(lagrange<1>(), BlockedInterleaved()));
  Dune::BlockVector<Ikarus::RealTuple<double, 2>> vBlocked(basis.size());

  auto localView = basis.localView();
  for (auto& ele :elements(gridView)) {
    localView.bind(ele);
    auto localBasis=Ikarus::LocalBasis(localView.tree().child(0).finiteElement().localBasis()) ;
    auto localF = Ikarus::SimpleLocalFunction(localBasis)
  }

}
