//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

//#include <fstream>
#include <vector>

#include <dune/typetree/leafnode.hh>
//#include <Eigen/Core>
#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/typetree/powernode.hh>

GTEST_TEST(Basis, Basistest) {
  auto gridFactory  = Dune::GridFactory<Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>>();
  const double L    = 1;
  const double h    = 1;
  const size_t elex = 2;
  const size_t eley = 1;
  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({L / 2, 0});
  gridFactory.insertVertex({L, 0});
  gridFactory.insertVertex({0, h});
  gridFactory.insertVertex({L / 2, h});
  gridFactory.insertVertex({L, h});

  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 3, 4});
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {1, 2, 4, 5});

  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Indices;
  constexpr int p = 1;
  auto basis      = makeBasis(gridView,
                         composite(power<2>(lagrange<p>(), FlatInterleaved()), lagrange<p - 1>(), FlatLexicographic()));

  auto dispBasis               = subspaceBasis(basis, _0);
  auto pressureBasis           = subspaceBasis(basis, _1);
  auto localView               = basis.localView();
  auto localViewOfDisplacement = dispBasis.localView();
  auto localViewOfPressure     = pressureBasis.localView();

  for (auto& e : elements(gridView)) {
    localView.bind(e);
    localViewOfDisplacement.bind(e);
    localViewOfPressure.bind(e);

    EXPECT_EQ(localView.size(), 9);                         // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localViewOfDisplacement.size(), 9);           // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localViewOfPressure.size(), 9);               // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localView.tree().size(), 9);                  // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localViewOfDisplacement.tree().size(), 8);    // Displacement Ansatzfunctions (DispDofs)
    EXPECT_EQ(localViewOfPressure.tree().size(), 1);        // Pressure Ansatzfunctions (PressureDofs)
    EXPECT_EQ(localViewOfDisplacement.tree().degree(), 2);  // How many displacement childs are there?
    EXPECT_EQ(localViewOfPressure.tree().degree(), 0);      // How many pressure childs are there?
    EXPECT_EQ(localViewOfDisplacement.tree().child(0).degree(), 0);
    EXPECT_EQ(localViewOfDisplacement.tree().child(1).degree(), 0);
    EXPECT_EQ(localViewOfDisplacement.globalBasis().dimension(), 14);  // How many dofs do we have in total?
    EXPECT_EQ(basis.dimension(), basis.size());                        // size and dimension cooincide for Flat basis
  }
}
