//
// Created by Alex on 21.07.2021.
//
#include <config.h>

#include <catch2/catch_test_macros.hpp>

#include "testHelpers.hh"

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <Eigen/Core>

#include <ikarus/solver/linearSolver/geometricMultigrid/gridTransfer.h>
#include <ikarus/utils/drawing/griddrawer.hh>
TEST_CASE("MultiGrid: TransferOperator", "[TransferOperatorTest.cpp]") {
  constexpr int gridDim = 2;
  using Grid            = Dune::YaspGrid<gridDim>;
  const double L        = 2;
  const double h        = 2;
  const size_t elex     = 1;
  const size_t eley     = 1;

  Dune::FieldVector<double, gridDim> bbox = {L, h};
  std::array<int, gridDim> eles           = {elex, eley};
  auto grid                               = std::make_shared<Grid>(bbox, eles);
  auto gridView                           = grid->leafGridView();

  grid->globalRefine(1);

  constexpr int numDofPerNode = 2;

  using namespace Dune::Functions::BasisFactory;
  auto preBasisFactory = power<gridDim>(lagrange<1>(), FlatInterleaved());

  Ikarus::GridTransfer transfer(grid);

  transfer.createOperators(preBasisFactory);

  Eigen::MatrixXd PExpected;
  PExpected.resize(numDofPerNode * grid->size(1, gridDim), numDofPerNode * grid->size(0, gridDim));
  PExpected << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0.25, 0,
      0.25, 0, 0.25, 0, 0.25, 0, 0, 0.25, 0, 0.25, 0, 0.25, 0, 0.25, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0,
      0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;

  Eigen::VectorXd v = Eigen::VectorXd::Random(numDofPerNode * grid->size(0, gridDim));
  Eigen::VectorXd vfine;
  transfer.prolongateFrom(0, v, vfine);

  CHECK_THAT(vfine, EigenApproxEqual(PExpected * v, 1e-15));
}

TEST_CASE("MultiGrid: TransferOperatorUGGridLshape", "[TransferOperatorTest.cpp]") {
  using namespace Ikarus;
  constexpr int gridDim = 2;
  //  //  /// ALUGrid Example
  using Grid = Dune::UGGrid<gridDim>;
  Dune::GridFactory<Grid> gridFactory;

  const double L = 1.0;

  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({L / 2, 0});
  gridFactory.insertVertex({0, L / 2});
  gridFactory.insertVertex({L / 2, L / 2});
  gridFactory.insertVertex({0, L});
  gridFactory.insertVertex({L / 2, L});
  gridFactory.insertVertex({L, L});
  gridFactory.insertVertex({L, L / 2});

  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 2, 3});
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {2, 3, 4, 5});
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {3, 7, 5, 6});

  auto grid = gridFactory.createGrid();
  grid->globalRefine(1);
  auto leafGridView = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;

  auto preBasisFactory = power<gridDim>(lagrange<1>(), FlatInterleaved());

  Ikarus::GridTransfer transfer(grid);

  transfer.createOperators(preBasisFactory);

  auto coarseBasis = makeBasis(grid->levelGridView(0), preBasisFactory);
  Eigen::VectorXd dCoarse(coarseBasis.size());
  for (int i = 0; i < dCoarse.size(); ++i) {
    dCoarse[i] = i;
  }

  Eigen::VectorXd dFine, dFineExpected(16);
  transfer.prolongateFrom(0, dCoarse, dFine);
  dFineExpected << 0, 1, 2, 3, 6, 7, 4, 5, 1, 2, 4, 5, 5, 6, 2, 3, 3, 4;

  CHECK_THAT(dFine, EigenApproxEqual(dFineExpected, 1e-15));

    Eigen::VectorXd dCoarse2;
    transfer.restrictTo(0,dFine,dCoarse2);

    CHECK_THAT(dCoarse, EigenApproxEqual(dCoarse2, 1e-15));
}